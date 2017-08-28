#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, subprocess
import numpy as np
import xmlschema
from ase.data import chemical_symbols, atomic_masses
from ase.dft.band_structure import BandStructure
from ase.calculators.calculator import all_changes, FileIOCalculator, Calculator, kpts2ndarray
import ase.units as units
from .io import get_atoms_from_xml_output


# Fix python3 types
try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str, bytes)


def get_band_structure(atoms=None, calc=None, ref=0):
    """
    Create band structure object from Atoms or calculator.

    This functions is rewritten here to allow the user to set the reference energy level.
    (The method calc.get_fermi_level() can fail in some cases as for insulators or for non-scf calculations)
    """

    atoms = atoms if atoms is not None else calc.atoms
    calc = calc if calc is not None else atoms.calc

    kpts = calc.get_bz_k_points()

    energies = []
    for s in range(calc.get_number_of_spins()):
        energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                         for k in range(len(kpts))])
    energies = np.array(energies)

    return BandStructure(cell=atoms.cell,
                         kpts=kpts,
                         energies=energies,
                         reference=ref)


class PostqeCalculator(Calculator):
    """
    This is a limited implementation of an ASE calculator for postqe.

    It does NOT generate an input file for espresso, NOR does it run pw.
    It reads some properties from an espresso xml output file (written
    according to espresso xml schema). The implemented properties are
    those needed for postprocessing.
    """

    implemented_properties = ['energy', 'forces']
    command = None

    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, command=None, schema=None, outdir=None, **kwargs):
        """
        File-IO calculator.

        :param restart:
        :param ignore_bad_restart_file:
        :param label:
        :param atoms:
        :param command: Command used to start calculation.
        :param schema:
        :param kwargs:
        """
        self.species = None
        self.schema = schema
        self.output = None
        self.outdir = outdir
        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, **kwargs)

    def calculate(self, atoms=None, properties=('energy',), system_changes=()):
        """
        This calculate method only reads the results from an xml output file
        in this calculator. The xml file is determined from the label field.

        :param atoms:
        :param properties:
        :param system_changes:
        :return:
        """
        raise NotImplementedError("Only for QE results post-processing")

    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def read_results(self):
        filename = self.label + '.xml'
        self.output = xmlschema.to_dict(filename, schema=self.schema, path=None)['output']
        self.atoms = get_atoms_from_xml_output(filename, output=self.output)
        self.results['energy'] = float(self.output["total_energy"]["etot"]) * units.Ry

    def band_structure(self, reference=0):
        """Create band-structure object for plotting.
        This method is redefined here to allow the user to set the reference energy (relying on the method get_fermi_level() is not safe.
        """
        return get_band_structure(calc=self, ref=reference)

    def get_number_of_bands(self):
        """Return the number of bands."""
        return int(self.output["band_structure"]["nbnd"])

    def get_xc_functional(self):
        """Return the XC-functional identifier.

        'LDA', 'PBE', ..."""
        return self.output["dft"]["functional"]

    def get_k_points(self):
        """Return all the k-points exactely as in the calculation.
        """
        nks = int(self.output["band_structure"]["nks"])  # get the number of k-points
        kpoints = np.zeros((nks, 3))
        ks_energies = self.output["band_structure"]["ks_energies"]

        for i in range(0, nks):
            for j in range(0, 3):
                kpoints[i, j] = float(ks_energies[i]['k_point']['$'][j])

        return kpoints

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2.
        """
        if self.output["magnetization"]["lsda"]:
            return 2
        return 1

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        if self.output["magnetization"]["lsda"]:
            return True
        return False

    def get_k_point_weights(self):
        """Weights of the k-points.

        The sum of all weights is one."""
        nks = int(self.output["band_structure"]["nks"])  # get the number of k-points
        weights = np.zeros(nks)
        ks_energies = self.output["band_structure"]["ks_energies"]

        for i in range(0, nks):
            weights[i] = float(ks_energies[i]['k_point']['@weight'])

        return weights

    def get_fermi_level(self):
        """Return the Fermi level.
        Warning:
        """
        try:
            ef = float(self.output["band_structure"]["fermi_energy"]) * units.Ry
            return ef
        except:
            raise Warning("Fermi energy not defined or not in output file")

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalues array."""

        nat = (self.output["atomic_structure"]["@nat"])
        nbnd = int(self.output["band_structure"]["nbnd"])
        ks_energies = (self.output["band_structure"]["ks_energies"])

        if self.get_spin_polarized():  # magnetic
            eigenvalues = np.zeros((nbnd//2))
            if spin == 0:
                # get bands for spin up
                for j in range(0, nbnd // 2):
                    # eigenvalue at k-point kpt, band j, spin up
                    eigenvalues[j] = float(ks_energies[kpt]['eigenvalues'][j]) * 2 * nat * units.Ry
            else:
                # get bands for spin down
                for j in range(nbnd // 2, nbnd):
                    # eigenvalue at k-point kpt, band j, spin down
                    eigenvalues[j - nbnd // 2] = float(ks_energies[kpt]['eigenvalues'][j]) * 2 * nat * units.Ry
        else:
            # non magnetic
            eigenvalues = np.zeros(nbnd)
            for j in range(0, nbnd):
                # eigenvalue at k-point kpt, band j
                eigenvalues[j] = float(ks_energies[kpt]['eigenvalues'][j]) * 2 * nat * units.Ry

        return eigenvalues

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""

        nbnd = int(self.output["band_structure"]["nbnd"])
        ks_energies = (self.output["band_structure"]["ks_energies"])

        if self.get_spin_polarized():
            # magnetic
            occupations = np.zeros((nbnd // 2))
            if spin == 0:
                # get bands for spin up
                for j in range(0, nbnd // 2):
                    # eigenvalue at k-point kpt, band j, spin up
                    occupations[j] = float(ks_energies[kpt]['occupations'][j])
            else:
                # get bands for spin down
                for j in range(nbnd // 2, nbnd):
                    # eigenvalue at k-point kpt, band j, spin down
                    occupations[j - nbnd // 2] = float(ks_energies[kpt]['occupations'][j])
        else:
            # non magnetic
            occupations = np.zeros(nbnd)
            for j in range(0, nbnd):
                # eigenvalue at k-point kpt, band j
                occupations[j] = float(ks_energies[kpt]['occupations'][j])

        return occupations

    # TODO: methods below are not implemented yet (do it if necessary)
    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal lattice vectors."""
        raise NotImplementedError
        #return kpoints

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal lattice vectors."""
        raise NotImplementedError
        #return kpoints

    def get_pseudo_density(self, spin=None, pad=True):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        raise NotImplementedError

    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        raise NotImplementedError

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        raise NotImplementedError


###############################################################################
###############################################################################
#
# In the following another calculator is implemented (as derived class),
#  which includes writing the input for QE and running the code with the old
#  text input format.
#
###############################################################################


# These lists contain the input parameters of pw.x divided per section (control,system,etc.)
# They correspond to Fortran namelists and are necessary to write the different sections of the input file from
# a unique dictionary of parameters (as in the base class FileIOCalculator)

control_keys = ['calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress',
                'tprnfor', 'dt', 'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr',
                'forc_conv_thr', 'disk_io', 'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm',
                'lberry', 'gdir', 'nppstr', 'lfcpopt', 'monopole']

system_keys = ['ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'tot_charge',
               'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3',
               'nr1s', 'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac',
               'occupations', 'one_atom_occupations', 'starting_spin_angle', 'degauss', 'smearing', 'nspin', 'noncolin',
               'ecfixed', 'qcutz', 'q2sigma', 'input_dft', 'exx_fraction', 'screening_parameter', 'exxdiv_treatment',
               'x_gamma_extrapolation', 'ecutvcut', 'nqx1', 'nqx2', 'nqx3', 'lda_plus_u', 'lda_plus_u_kind',
               'Hubbard_U', 'Hubbard_J0', 'Hubbard_alpha', 'Hubbard_beta', 'Hubbard_J(i,ityp)',
               'starting_ns_eigenvalue(m,ispin,I)', 'U_projection_type', 'edir', 'emaxpos', 'eopreg', 'eamp',
               'angle1', 'angle2', 'constrained_magnetization', 'fixed_magnetization', 'lambda', 'report', 'lspinorb',
               'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit', 'fcp_mu', 'vdw_corr',
               'london', 'london_s6', 'london_c6', 'london_rvdw', 'london_rcut', 'ts_vdw_econv_thr', 'ts_vdw_isolated',
               'xdm', 'xdm_a1', 'xdm_a2', 'space_group', 'uniqueb', 'origin_choice', 'rhombohedral', 'zmon', 'realxz',
               'block', 'block_1', 'block_2', 'block_height']

electrons_keys = ['electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init',
                  'conv_thr_multi', 'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns', 'diagonalization',
                  'ortho_para', 'diago_thr_init', 'diago_cg_maxiter', 'diago_david_ndim', 'diago_full_acc', 'efield',
                  'efield_cart', 'efield_phase', 'startingpot', 'startingwfc', 'tqr']


ions_keys = ['ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot',
             'ion_temperature', 'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim',
             'trust_radius_max', 'trust_radius_min', 'trust_radius_ini', 'w_1', 'w_2']

cell_keys = ['cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree']


def write_type(f, key, value):
    """ Write a line "key = value,\n" in file f according to value type"""
    if isinstance(value, bool):
        f.write('    %s = .%s.,\n' % (key, str(value).lower()))
    elif isinstance(value, float):
        f.write('    %s = %g,\n' % (key, value))
    elif isinstance(value, int):
        f.write('    %s = %d,\n' % (key, value))
    elif isinstance(value, basestring):
        f.write("    %s = '%s',\n" % (key, value))


class EspressoCalculator(PostqeCalculator):
    """
    This is a full calculator for Quantum Espresso, including generating input file and running the code.
    The input file is the old text format which will be substituted by a new xml format in the future.

    To use the full calculator, some input parameters for QE must be handled.
    """
    implemented_properties = ['energy', 'forces']
    command = '/home/brunato/Development/projects/postqe/postqe/fortran/build/q-e/bin/pw.x < PREFIX.in > PREFIX.out'

    # These are reasonable default values for the REQUIRED parameters in pw.x input.
    # All other default values are not set here and let to pw.x, unless the user defines them in the calculator.
    default_parameters = dict(
        ibrav=0,
        nat=1,
        ntyp=1,
        ecutwfc=50,
        kpoints=[1, 1, 1, 0, 0, 0]
    )

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='pw', atoms=None, command=None, pp_dict=None, **kwargs):
        """Construct Espresso-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.out, ...).
            Default is 'pw'.

        Examples
        ========
        Use default values:

        >>> from ase import Atoms
        >>> h = Atoms('H', calculator=EspressoCalculator(ecut=200, toldfe=0.001))
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()

        """
        # self.species = None
        self.pp_dict = pp_dict

        PostqeCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)

        if command is not None:
            self.command = command
        else:
            name = 'ASE_' + self.name.upper() + '_COMMAND'
            self.command = os.environ.get(name, self.command)

    def check_state(self, atoms, tol=1e-15):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input parameters to files-file."""
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        if 'numbers' in system_changes or 'initial_magmoms' in system_changes:
            self.initialize(atoms)

        print ("Writing input file... " + self.label + ".in")
        param = self.parameters   # copy the parameters into param

        finput = open(self.label + '.in', 'w')

        # Write CONTROL section
        finput.write(' &CONTROL \n')
        for tag in param:
            if tag in control_keys:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write SYSTEM section
        finput.write(' &SYSTEM \n')
        for tag in param:
            if tag in system_keys:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write ELECTRONS section
        finput.write(' &ELECTRONS \n')
        for tag in param:
            if tag in electrons_keys:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write IONS section
        finput.write(' &IONS \n')
        for tag in param:
            if tag in ions_keys:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write CELL section
        finput.write(' &CELL \n')
        for tag in param:
            if tag in ions_keys:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write the ATOMIC_SPECIES section
        finput.write('ATOMIC_SPECIES\n')
        for Z in self.species:
            finput.write(chemical_symbols[Z] + ' ' + str(atomic_masses[Z]) + ' ')
            finput.write('%s\n' % self.pp_dict[chemical_symbols[Z]])
            pass

        # Write the ATOMIC_POSITIONS section
        # TODO check positions are in Bohr or Ang or else
        finput.write('ATOMIC_POSITIONS {Bohr}\n')
        for i, pos in zip(atoms.numbers, atoms.positions):
            finput.write('%s ' % chemical_symbols[i])
            finput.write('%.14f %.14f %.14f\n' % tuple(pos))

        # Write the CELL_PARAMETERS section (assume they are in Angstrom)
        finput.write('CELL_PARAMETERS {Ang}\n')
        for v in atoms.cell:
            finput.write('%.14f %.14f %.14f\n' % tuple(v))

        # Write the K_POINT section, always as a list of k-points
        # For Monkhorst meshes, the method "calculate" has already generated a list of k-points in
        # self.kpts using ASE function kpts2ndarray
        finput.write('K_POINTS tpiba\n')
        finput.write('%d\n' % len(self.kpts))   # first write the number of k-points
        for i in range(0, len(self.kpts)):
            finput.write('%f %f %f' % tuple(self.kpts[i]) + ' 1.0\n')   # assume unary weight for all

        # TODO: check if it makes sense in some cases to let QE generate the Monkhorst mesh
        # finput.write('K_POINTS AUTOMATIC\n')
        # finput.write('%d %d %d %d %d %d\n' % tuple(param['kpoints']))

        finput.close()

    def initialize(self, atoms):

        self.parameters['nat'] = atoms.get_number_of_atoms()  # how many atoms are in the atoms object
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)
        self.parameters['ntyp'] = len(self.species)
        self.spinpol = atoms.get_initial_magnetic_moments().any()

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.kpts = kpts2ndarray(self.parameters.kpts, atoms)
        self.write_input(self.atoms, properties, system_changes)

        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('ASE_' + self.name.upper() + '_COMMAND') +
                               'or supply the command keyword')
        command = self.command.replace('PREFIX', self.prefix)
        olddir = os.getcwd()
        try:
            print(self.directory)
            os.chdir(self.directory)
            print (command)
            errorcode = subprocess.call(command, shell=True)
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s in %s returned an error: %d' %
                               (self.name, self.directory, errorcode))

        self.read_results()

    def read_results(self):
        # FIXME: compose XML output filename from parameters
        filename = self.directory + '/temp/pwscf.xml'
        self.output = xmlschema.to_dict(filename, schema=self.schema, path="./output")
        self.atoms = get_atoms_from_xml_output(filename, output=self.output)
        self.results['energy'] = float(self.output["total_energy"]["etot"]) * units.Ry
