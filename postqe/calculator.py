#
# Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#
from distutils.log import error
from importlib.resources import path
import os
import re
import subprocess
from tabnanny import check
import numpy as np
import pathlib
import qeschema

from ase.atoms import Atoms, Atom
from ase.data import chemical_symbols, atomic_masses
from ase.spectrum.band_structure import BandStructure
from ase.calculators.calculator import all_changes, FileIOCalculator, Calculator, kpts2ndarray
import ase.units as units


# These tuples contain the input parameters of pw.x divided per section (control, system, etc.)
# They correspond to Fortran namelists and are necessary to write the different sections of
# the input file from a unique dictionary of parameters (as in the base class FileIOCalculator)

CONTROL_KEYS = ('calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect',
                'nstep', 'iprint', 'tstress', 'tprnfor', 'dt', 'outdir', 'wfcdir',
                'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr',
                'disk_io', 'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc',
                'lorbm', 'lberry', 'gdir', 'nppstr', 'lfcpopt', 'monopole')

SYSTEM_KEYS = ('ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat',
               'ntyp', 'nbnd', 'tot_charge', 'tot_magnetization', 'starting_magnetization',
               'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3', 'nr1s', 'nr2s', 'nr3s',
               'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac',
               'occupations', 'one_atom_occupations', 'starting_spin_angle', 'degauss',
               'smearing', 'nspin', 'noncolin', 'ecfixed', 'qcutz', 'q2sigma', 'input_dft',
               'exx_fraction', 'screening_parameter', 'exxdiv_treatment',
               'x_gamma_extrapolation', 'ecutvcut', 'nqx1', 'nqx2', 'nqx3', 'lda_plus_u',
               'lda_plus_u_kind', 'Hubbard_U', 'Hubbard_J0', 'Hubbard_alpha', 'Hubbard_beta',
               'Hubbard_J(i,ityp)', 'starting_ns_eigenvalue(m,ispin,I)', 'U_projection_type',
               'edir', 'emaxpos', 'eopreg', 'eamp', 'angle1', 'angle2',
               'constrained_magnetization', 'fixed_magnetization', 'lambda', 'report',
               'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit',
               'fcp_mu', 'vdw_corr', 'london', 'london_s6', 'london_c6', 'london_rvdw',
               'london_rcut', 'ts_vdw_econv_thr', 'ts_vdw_isolated', 'xdm', 'xdm_a1',
               'xdm_a2', 'space_group', 'uniqueb', 'origin_choice', 'rhombohedral',
               'zmon', 'realxz', 'block', 'block_1', 'block_2', 'block_height')

ELECTRONS_KEYS = ('electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr',
                  'conv_thr_init', 'conv_thr_multi', 'mixing_mode', 'mixing_beta',
                  'mixing_ndim', 'mixing_fixed_ns', 'diagonalization', 'ortho_para',
                  'diago_thr_init', 'diago_cg_maxiter', 'diago_david_ndim',
                  'diago_full_acc', 'efield', 'efield_cart', 'efield_phase',
                  'startingpot', 'startingwfc', 'tqr')


IONS_KEYS = ('ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation',
             'remove_rigid_rot', 'ion_temperature', 'tempw', 'tolp', 'delta_t', 'nraise',
             'refold_pos', 'upscale', 'bfgs_ndim', 'trust_radius_max', 'trust_radius_min',
             'trust_radius_ini', 'w_1', 'w_2')

CELL_KEYS = ('cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree')


def write_type(f, key, value):
    """ Write a line "key = value,\n" in file f according to value type"""
    if isinstance(value, bool):
        f.write('    %s = .%s.,\n' % (key, str(value).lower()))
    elif isinstance(value, float):
        f.write('    %s = %g,\n' % (key, value))
    elif isinstance(value, int):
        f.write('    %s = %d,\n' % (key, value))
    elif isinstance(value, str):
        f.write("    %s = '%s',\n" % (key, value))


def split_atomic_symbol(x):
    regex = r"([a-zA-Z]{1,2})([1-9]?\d?)"
    match = re.match(regex, x)
    if match:
        return match.groups()
    else:
        return False


def get_band_structure(atoms=None, calc=None, ref=0):
    """
    Create band structure object from Atoms or calculator.

    This functions is rewritten here to allow the user to set the
    reference energy level. The method calc.get_fermi_level() can
    fail in some cases as for insulators or for non-scf calculations.
    """

    atoms = atoms if atoms is not None else calc.atoms
    calc = calc if calc is not None else atoms.calc

    kpts = calc.get_k_points()

    energies = []
    for s in range(calc.get_number_of_spins()):
        energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                         for k in range(len(kpts))])
    energies = np.array(energies)

    return BandStructure(path=atoms.cell.bandpath(path=kpts),
                         energies=energies,
                         reference=ref)


class EspressoCalculator(FileIOCalculator):
    """
    An ASE calculator for Quantum Espresso, including generating input file
    and running the code. The input file for the current release is the old
    text format which will be substituted by a new XML format in the future.

    :param restart: prefix for restart file. May contain a directory. Default \
    is `None`: don't restart.
    :param ignore_bad_restart_file: ignore broken or missing restart file. \
    For default it is an error if the restart file is missing or broken.
    :param label: Prefix used for all files, as like as PW's 'prefix' parameter. \
    Defaults to 'pwscf'.
    :param atoms: optional Atoms object to which the calculator will be attached. \
    When restarting, atoms will get its positions and unit-cell updated from file.
    :param command: Command used to start calculation.
    :param outdir: directory containing the input data. Default to the value of \
    ESPRESSO_TMPDIR environment variable if set or current directory ('.') otherwise
    :param schema:
    :param kwargs:

    Examples
    ========
    Use default values:

    >>> from ase import Atoms
    >>> h = Atoms('H', calculator=EspressoCalculator(ecut=200, toldfe=0.001))
    >>> h.center(vacuum=3.0)
    >>> e = h.get_potential_energy()
    """
    ignored_changes = {'pbc'}
    discard_results_on_any_change = True

    implemented_properties = ['energy', 'forces']

    command = str(pathlib.Path(__file__).parent.joinpath(
        'fortran/build/q-e/bin/pw.x < PREFIX.in > PREFIX.out'
    ))

    # These are reasonable default values for the REQUIRED parameters in pw.x input.
    # All other default values are not set here and let to pw.x, unless the user
    # defines them in the calculator.
    default_parameters = dict(
        ibrav=0,
        nat=1,
        ntyp=1,
        ecutwfc=50,
        kpoints=[1, 1, 1, 0, 0, 0]
    )
    species = None
    spinpol = None
    kpts = None

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='pwscf', atoms=None, command=None, outdir=None,
                 schema=None, pp_dict=None, **kwargs):

        # Check outdir and label to fit the preferred format for an ASE calculator
        if outdir is not None:
            pass
        elif '/' in label:
            outdir, label = label.rsplit('/', 1)
        else:
            outdir = os.environ.get('ESPRESSO_TMPDIR', '.')

        if not os.path.isdir(outdir):
            raise ValueError("{!r} is not a directory path".format(outdir))

        self.pp_dict = pp_dict
        self.xml_document = qeschema.PwDocument(schema=schema)
        super().__init__(restart, ignore_bad_restart_file, label, atoms,
                         command, directory=outdir, **kwargs)

    @property
    def outdir(self):
        return self.directory

    @property
    def schema(self):
        return self.xml_document.schema

    @property
    def input(self):
        return self.xml_document.to_dict(preserve_root=False)['input']

    @property
    def output(self):
        return self.xml_document.to_dict(preserve_root=False)['output']

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input parameters to files-file."""
        super().write_input(atoms, properties, system_changes)

        if 'numbers' in system_changes or 'initial_magmoms' in system_changes:
            self.initialize(atoms)

        print("Writing input file... " + self.label + ".in")
        param = self.parameters   # copy the parameters into param

        finput = open(self.label + '.in', 'w')

        # Write CONTROL section
        finput.write(' &CONTROL \n')
        for tag in param:
            if tag in CONTROL_KEYS:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write SYSTEM section
        finput.write(' &SYSTEM \n')
        for tag in param:
            if tag in SYSTEM_KEYS:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write ELECTRONS section
        finput.write(' &ELECTRONS \n')
        for tag in param:
            if tag in ELECTRONS_KEYS:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write IONS section
        finput.write(' &IONS \n')
        for tag in param:
            if tag in IONS_KEYS:
                write_type(finput, tag, param[tag])
        finput.write('/\n')

        # Write CELL section
        finput.write(' &CELL \n')
        for tag in param:
            if tag in IONS_KEYS:
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
            finput.write('%f %f %f' % tuple(self.kpts[i]) + ' 1.0\n')  # assume unary weight for all

        # TODO: check if it makes sense in some cases to let QE generate the Monkhorst mesh
        # finput.write('K_POINTS AUTOMATIC\n')
        # finput.write('%d %d %d %d %d %d\n' % tuple(param['kpoints']))

        finput.close()

    def initialize(self, atoms):
        self.parameters['nat'] = atoms.get_number_of_atoms()
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)
        self.parameters['ntyp'] = len(self.species)
        self.spinpol = atoms.get_initial_magnetic_moments().any()

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
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
            print(command)
            errorcode = subprocess.call(command, shell=True)
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s in %s returned an error: %d' %
                               (self.name, self.directory, errorcode))

        self.read_results()

    def read_results(self, filename=None):
        if filename is None:
            filename = str(pathlib.Path(self.label).joinpath('data-file-schema.xml'))
        self.xml_document.read(filename)
        self.atoms = self.get_atoms_from_xml_output()
        self.results['energy'] = float(self.output["total_energy"]["etot"]) * units.Ry

    def band_structure(self, reference=0):
        """Create band-structure object for plotting.

        This method is redefined here to allow the user to set the reference
        energy (relying on the method get_fermi_level() is not safe).
        """
        return get_band_structure(calc=self, ref=reference)

    def get_number_of_bands(self):
        """Return the number of bands."""
        return int(self.output["band_structure"]["nbnd"])

    def get_xc_functional(self):
        """Returns the XC-functional identifier ('LDA', 'PBE', ...)."""
        return self.output["dft"]["functional"]

    def get_k_points(self):
        """Returns all the k-points exactly as in the calculation."""
        nks = int(self.output["band_structure"]["nks"])  # get the number of k-points
        kpoints = np.zeros((nks, 3))
        ks_energies = self.output["band_structure"]["ks_energies"]

        for i in range(0, nks):
            for j in range(0, 3):
                kpoints[i, j] = float(ks_energies[i]['k_point']['$'][j])

        return kpoints

    def get_atoms_from_xml_output(self):
        """
        Returns an Atoms object constructed from an XML QE file (according to the schema).

        :return: An Atoms object.
        """
        atomic_structure = self.output['atomic_structure']
        a1 = np.array(atomic_structure["cell"]["a1"])
        a2 = np.array(atomic_structure["cell"]["a2"])
        a3 = np.array(atomic_structure["cell"]["a3"])
        a_p = (atomic_structure["atomic_positions"]["atom"])

        atoms = Atoms()

        # First define the unit cell from a1, a2, a3 and alat
        cell = np.zeros((3, 3))
        cell[0] = a1
        cell[1] = a2
        cell[2] = a3
        atoms.set_cell(cell)

        # Now the atoms in the unit cell
        for atomx in a_p:
            # TODO: extend to all possible cases the symbol splitting
            #  (for now, only numbering up to 9 work). Not a very common case...
            symbol = split_atomic_symbol(atomx['@name'])[0]
            x = float(atomx['$'][0])
            y = float(atomx['$'][1])
            z = float(atomx['$'][2])
            atoms.append(Atom(symbol, (x, y, z)))
            atoms.pbc=[1,1,1]

        return atoms

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2.
        """
        return 1 if not self.get_spin_polarized() else 2

    def get_spin_polarized(self):
        """Test if it is a spin-polarized calculation."""
        return self.output["magnetization"]["lsda"] is True

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
            ef = float(self.output["band_structure"]["fermi_energy"]) * units.Ha
            return ef
        except (TypeError, KeyError):
            raise Warning("Fermi energy not defined or not in output file")

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalues array.
           For spin polarized specify spin=1 or spin=2, default spin=1
        """

        nat = (self.output["atomic_structure"]["@nat"])
        try:
            nbnd = int(self.output["band_structure"]["nbnd"])
            use_updw = False
        except KeyError:
            nbnd_up = int(self.output["band_structure"]["nbnd_up"])
            nbnd_dw = int(self.output["band_structure"]["nbnd_dw"])
            use_updw = True 
        ks_energies = (self.output["band_structure"]["ks_energies"])

        if self.get_spin_polarized():  # magnetic
            if spin == 0:
                spin = 1 
            if not use_updw:
                nbnd_up = nbnd // 2 
                nbnd_dw = nbnd // 2 
            if spin == 1:
                # get bands for spin up 
                eigenvalues = np.zeros(nbnd_up)
                for j in range(0, nbnd_up):
                    # eigenvalue at k-point kpt, band j, spin up
                    try:
                        eigenvalues[j] = float(ks_energies[kpt]['eigenvalues'][j])  * units.Ha
                    except KeyError:
                        eigenvalues[j] = float(ks_energies[kpt]['eigenvalues']['$'][j]) * units.Ha
            else:
                # get bands for spin down
                eigenvalues = np.zeros(nbnd_dw)
                for j in range(nbnd_up, nbnd_up + nbnd_dw ):
                    # eigenvalue at k-point kpt, band j, spin down
                    try:
                        eigenvalues[j - nbnd_up] = float(ks_energies[kpt]['eigenvalues'][j]) * units.Ha
                    except KeyError:
                        eigenvalues[j - nbnd_up] = float(ks_energies[kpt]['eigenvalues']['$'][j]) * units.Ha
        else:
            # non magnetic
            eigenvalues = np.zeros(nbnd)
            for j in range(0, nbnd):
                # eigenvalue at k-point kpt, band j
                try:
                    eigenvalues[j] = float(ks_energies[kpt]['eigenvalues'][j]) * units.Ha
                except KeyError:
                    eigenvalues[j] = float(ks_energies[kpt]['eigenvalues']['$'] [j]) * units.Ha

        return eigenvalues

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array.  For spin polarized case specify spin=1 or spin=2"""
        try:
            nbnd = int(self.output["band_structure"]["nbnd"])
            use_updw = False
        except KeyError:
            nbnd_up = int(self.output["band_structure"]["nbnd_up"])
            nbnd_dw = int(self.output["band_structure"]["nbnd_dw"])
            use_updw = True
        #
        ks_energies = (self.output["band_structure"]["ks_energies"])

        if self.get_spin_polarized():
            # magnetic
            if spin == 0:
                spin = 1 
            elif spin > 2: 
                raise ValueError("Spin can be either 1 or 2 ")
            if not use_updw:
                nbnd_up = nbnd // 2 
                nbnd_dw = nbnd // 2 
            if spin == 1:
                # get bands for spin up
                occupations = np.zeros(nbnd_up) 
                for j in range(0, nbnd_up):
                    # eigenvalue at k-point kpt, band j, spin up
                    try:
                        occupations[j] = float(ks_energies[kpt]['occupations'][j])
                    except KeyError:
                        occupations[j] = float(ks_energies[kpt]['occupations']['$'][j])  
            else:
                # get bands for spin down
                occupations = np.zeros(nbnd_dw)
                for j in range(nbnd_up, nbnd_up + nbnd_dw):
                    # eigenvalue at k-point kpt, band j, spin down
                    try:
                        occupations[j - nbnd_up ] = float(ks_energies[kpt]['occupations'][j])
                    except KeyError:
                        occupations[j - nbnd_up ] = float(ks_energies[kpt]['occupations']['$'][j])
        else:
            # non magnetic
            occupations = np.zeros(nbnd)
            for j in range(0, nbnd):
                # eigenvalue at k-point kpt, band j
                try: 
                    occupations[j] = float(ks_energies[kpt]['occupations'][j])
                except KeyError:
                    occupations[j] = float(ks_energies[kpt]['occupations']['$'][j])  

        return occupations

    # Here are a number of additional getter methods, some very specific for Quantum Espresso
    def get_nr(self):
        nr = np.array([self.output["basis_set"]["fft_grid"]["@nr1"],
                       self.output["basis_set"]["fft_grid"]["@nr2"],
                       self.output["basis_set"]["fft_grid"]["@nr3"]], int)
        return nr

    def get_nr_smooth(self):
        nr_smooth = np.array([self.output["basis_set"]["fft_smooth"]["@nr1"],
                              self.output["basis_set"]["fft_smooth"]["@nr2"],
                              self.output["basis_set"]["fft_smooth"]["@nr3"]], int)
        return nr_smooth

    def get_ibrav(self):
        return int(self.output["atomic_structure"]["@bravais_index"])

    def get_alat(self):
        return float(self.output["atomic_structure"]["@alat"])

    def get_ecutwfc(self):
        return float(self.output["basis_set"]["ecutwfc"])

    def get_ecutrho(self):
        return float(self.output["basis_set"]["ecutrho"])

    def get_pseudodir(self):
        return self.input["control_variables"]["pseudo_dir"]

    # TODO: these two methods are just a temporary patch
    #  (a and b vectors can be obtained from Atoms object)
    def get_a_vectors(self):
        a1 = np.array(self.output["atomic_structure"]["cell"]["a1"])
        a2 = np.array(self.output["atomic_structure"]["cell"]["a2"])
        a3 = np.array(self.output["atomic_structure"]["cell"]["a3"])
        return np.array([a1, a2, a3])

    def get_b_vectors(self):
        b1 = np.array(self.output["basis_set"]["reciprocal_lattice"]["b1"])
        b2 = np.array(self.output["basis_set"]["reciprocal_lattice"]["b2"])
        b3 = np.array(self.output["basis_set"]["reciprocal_lattice"]["b3"])
        return np.array([b1, b2, b3])

    def get_atomic_positions(self):
        atomic_positions = self.output["atomic_structure"]["atomic_positions"]["atom"]
        if not isinstance(atomic_positions, list):
            return [atomic_positions]
        return atomic_positions

    def get_atomic_species(self):
        atomic_species = self.output["atomic_species"]["species"]

        # for subsequent loops it is important to return a list.
        if not isinstance(atomic_species, list):
            return [atomic_species]
        return atomic_species

    # TODO: methods below are not implemented yet (do it if necessary)
    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.
        The coordinates are relative to reciprocal lattice vectors."""
        kpoints= self.get_k_points() / self.get_alat()   
        m = self.get_a_vectors() 
        res = kpoints.dot(m.T) 
        nint = lambda d: int(round(d,0))
        center = lambda x: round(x - nint(x), 8) 
        res = np.array( [np.array( [center(c) for c in _]) for _ in res[:]])
        return res

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal lattice vectors."""
        raise NotImplementedError
        # return kpoints

    def get_pseudo_density(self, dataset='total', direction=3, check_noncolin = False, filename=None):
        """Return pseudo-density array. 
           dataset string is used  to chose which density is retrieved, possible values are 
           'total' (default) and 'magnetization' For the noncollinear case the variable direction 
           specifies the direction (3 is the default)
        """
        from .charge import get_charge_r, get_magnetization_r 
        if filename is None:
            filename = str(pathlib.Path( self.label ).joinpath( 'charge-density.hdf5' ))
        if dataset == 'total':
            return get_charge_r(filename)
        elif dataset == 'magnetization':
            temp = get_magnetization_r (filename, direction = direction )
            if temp is None:
                return None 
            if check_noncolin:
                if not temp[0]:
                    raise Exception("Non collinear magnetization not found") 
            return temp[1]



    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        raise NotImplementedError

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, filename=None):
        """Return pseudo-wave-function array.  
        for data read from the default  directory  kpt and spin are needed. 
        for data read from file specified with filename, kpt and spin are neglected. 
        :band: integer specify what band has to be extracted 
        :kpt:  integer specify which k point is read. 
        :spin: 1 or 2 for the lsda case specifies  which spin has to be selected 
        :filename: string with the path to a custom wavefunction hdf5 file. """
        from .readutils import get_wf_attributes, get_wavefunctions, get_wfc_miller_indices 
        from .charge    import charge_r_from_cdata
        spinlabels = ['up','dw']
        if filename is not None:
            attrs  = get_wf_attributes(filename=filename)
            kpt = attrs['ik']
        else: 
            if self.get_spin_polarized():
                filename = str( pathlib.Path( self.label ).joinpath( "".join(['wfc', spinlabels[spin-1], str(kpt),'.hdf5'] )))
            else:
                filename =  str( pathlib.Path( self.label ).joinpath( "".join(['wfc', str(kpt),'.hdf5'] )))
            attrs = get_wf_attributes(filename = filename)
        #
        xk = attrs['xk']
        igwx = attrs['igwx']
        data = get_wavefunctions(filename,band-1,band)[0] 
        MI = get_wfc_miller_indices(filename)
        nr1 = 2*max(abs(MI[:, 0]))+1
        nr2 = 2*max(abs(MI[:, 1]))+1
        nr3 = 2*max(abs(MI[:, 2]))+1
        gamma_only = 'TRUE' in str(attrs['gamma_only']).upper() 
        nr = np.array([nr1, nr2, nr3])
        return charge_r_from_cdata(data, MI, gamma_only, nr) 



    

# noinspection PyAbstractClass
class PostqeCalculator(EspressoCalculator):
    """
    This is a restricted ASE file-IO calculator for Quantum Espresso.
    It does NOT generate an input file for Espresso, NOR does it run pw.
    It reads some properties from an Espresso XML output file (written
    according to espresso XSD schema). The implemented properties are
    those needed for postprocessing.
    """
    def __init__(self, restart=None, ignore_bad_restart_file=False, label='pwscf',
                 atoms=None, outdir=None, schema=None, **kwargs):
        kwargs.pop('command', None)
        kwargs.pop('pp_dict', None)
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file,
                         label=label, atoms=atoms, outdir=outdir, schema=schema, **kwargs)
        self.command = None

    def calculate(self, atoms=None, properties=('energy',), system_changes=()):
        raise RuntimeError("{} is only for QE results post-processing".format(type(self)))
