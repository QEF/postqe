#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import re
import platform
from pathlib import Path

from setuptools import setup
from setuptools.command.install_lib import install_lib

from distutils import log
from distutils.command.clean import clean  # type: ignore[attr-defined]
from distutils.file_util import move_file


VERSION_NUMBER_PATTERN = re.compile(r"version_number\s*=\s*(\'[^\']*\'|\"[^\"]*\")")

with Path(__file__).parent.joinpath("requirements.txt").open() as fp:
    REQUIREMENTS = fp.read()


def find_extension_module(pattern):
    """
    Returns the absolute pathname of an extension module built for the running platform.

    :param pattern: a glob pattern iterable for finding the candidate extension modules.
    :return: A pathname string or `None` if no suitable module is found.
    """
    python_version = "".join(platform.python_version_tuple()[:2])

    for filename in map(str, pattern):
        if (
            f"-{python_version}-" not in filename
            and f"{python_version}m" not in filename
        ):
            continue
        elif platform.python_implementation().lower() not in filename:
            continue
        elif platform.system().lower() not in filename:
            continue
        elif platform.machine() not in filename:
            continue
        return filename
    else:
        return None


class InstallLibCommand(install_lib):
    verbose: bool
    dry_run: bool

    def build(self):
        # base method builds python modules (build_py) and C extensions (build_ext)
        install_lib.build(self)

        if not self.skip_build:
            # Build Fortran based modules
            # TODO: clean steps, checking if binary modules already exist
            #   and reading the QE_TOPDIR variable
            build_dir = Path(self.build_dir).joinpath('postqe/fortran')
            qe_build_dir = build_dir.joinpath("qe")
            if not qe_build_dir.exists():
                qe_build_dir.mkdir(parents=True)

            self._build_f90utils_module(build_dir)
            self._build_quantum_espresso(qe_build_dir)
            self._build_pyqe_module(build_dir)

            move_file(
                src=str(build_dir.joinpath("pyqe.py")),
                dst=str(build_dir.parent),
                verbose=self.verbose,
                dry_run=self.dry_run,
            )

            for srcfile in build_dir.glob("_pyqe*.so"):
                move_file(
                    src=str(srcfile),
                    dst=str(build_dir.parent),
                    verbose=self.verbose,
                    dry_run=self.dry_run,
                )

            for srcfile in build_dir.glob("f90utils*.so"):
                move_file(
                    src=str(srcfile),
                    dst=str(build_dir.parent),
                    verbose=self.verbose,
                    dry_run=self.dry_run,
                )

    @staticmethod
    def _build_f90utils_module(build_dir):
        print("Build f90utils module ...")
        os.system(f"make -f f90utils_Makefile -C {build_dir} f90utils_module")

    @staticmethod
    def _build_quantum_espresso(qe_build_dir):
        try:
            qe_topdir = Path(os.environ["QE_TOPDIR"]).absolute()
        except KeyError:
            print("Please specify the QE_TOPDIR before executing setup.py")
            raise
        else:
            if not qe_topdir.is_dir():
                raise FileNotFoundError(
                    f"Missing QE_TOPDIR={str(qe_topdir)!r} directory!"
                )

        # Check QE installation
        for version_file in qe_topdir.glob("include/*version.h"):
            with version_file.open() as _fp:
                version_number = VERSION_NUMBER_PATTERN.search(_fp.read())

            if version_number is None:
                raise ValueError(f"Missing version label in {version_file}")

            print("QE {}".format(version_number.group(0)))
            if version_number.group(1).strip("'\"") < "7.1":
                raise ValueError("Build aborted: QE version too old!")

            break
        else:
            raise FileNotFoundError("Missing QE version file info!")

        print("Configure Quantum Espresso ...")
        os.system(
            f"cmake -DQE_ENABLE_MPI=off "
            f"-DCMAKE_Fortran_FLAGS=-fPIC "
            f"-DCMAKE_C_FLAGS=-fPIC "
            f"-B {qe_build_dir} "
            f"-S {os.environ['QE_TOPDIR']}"
        )

        print("Build Quantum Espresso ...")
        os.system(f"make -C {qe_build_dir} -j qe_pp")

    @staticmethod
    def _build_pyqe_module(build_dir):
        print("Build pyqe module ...")

        qe_build_dir = build_dir.joinpath("qe").absolute()
        # Read the CMakeCache.txt file to find:
        # LAPACK LIBS path
        # BLAS LIBS path
        # FFTW3 LIBS path
        with open(f"{qe_build_dir}/CMakeCache.txt", "r") as file:
            print(file)
            cmake_cache_lines = list(filter(
                lambda s: 'FIND_PACKAGE_MESSAGE_DETAILS' in s,
                file.readlines()))

        def get_path(lib, lines):
            pattern = r"(\[[^\s\[\]]+\])(\[[^\s\[\]]+\])+"
            pattern2 = r"([/-][^;\s\[\]]+)"
            line = next(filter(lambda s: lib in s, lines), None)
            libpath = None if line is None else " ".join(
                re.findall(pattern2, re.search(pattern, line).group(1))
            )
            return libpath

        lapack_path = get_path('LAPACK', cmake_cache_lines)
        blas_path = get_path('BLAS', cmake_cache_lines)
        fftw3_first_path = get_path('FFTW', cmake_cache_lines)

        if fftw3_first_path:
            print("FFTW3 first path found:", fftw3_first_path)
            os.environ["FFTW3_LIBS"] = fftw3_first_path
        else:
            print("No FFTW3 path found")

        if lapack_path:
            print("LAPACK path found:", lapack_path)
            os.environ["LAPACK_LIBS"] = lapack_path
        else:
            print("No LAPACK path found")

        if blas_path:
            print("BLAS path found:", blas_path)
            if blas_path != lapack_path:
                os.environ["BLAS_LIBS"] = blas_path
        else:
            print("No BLAS path found")

        # Run the make target to build the pyqe module
        os.system(
            f"make LAPACK_LIBS=\"{os.environ.get('LAPACK_LIBS', '')}\" "
            f"BLAS_LIBS=\"{os.environ.get('BLAS_LIBS', '')}\" "
            f"FFTW3_LIBS=\"{os.environ.get('FFTW3_LIBS', '')}\" "
            f"QE_BUILD_DIR={qe_build_dir} "
            f"BUILD_DIR={build_dir.absolute()} "
            f"-f pyqe_Makefile "
            f"-C {build_dir} "
            f"pyqe"
        )

        # Modify python wrapper module, fixing import in postqe package.
        with build_dir.joinpath("pyqe.py").open() as _fp:
            python_wrapper_lines = _fp.readlines()

        with build_dir.joinpath("pyqe.py").open(mode="w") as _fp:
            python_wrapper_lines[0] = "# Altered wrapper for postqe\n"
            python_wrapper_lines[1] = "from . import _pyqe\n"
            _fp.writelines(python_wrapper_lines)


class CleanCommand(clean):
    def run(self):
        clean.run(self)
        if self.all:
            log.info("Remove built extension modules ...")
            source_dir = Path(__file__).parent.joinpath("postqe")

            wrapper_module = source_dir.joinpath("pyqe.py")
            if wrapper_module.is_file():
                log.info(f"Remove {wrapper_module}")
                wrapper_module.unlink()

            for ext_module in source_dir.glob("*.so"):
                log.info(f"Remove {ext_module}")
                os.unlink(ext_module)


setup(
    name='postqe',
    version='1.0.0',
    packages=['postqe'],
    package_data={
        'postqe': ['fortran/*', 'fortran/wrapfiles/*.f90']
    },
    install_requires=REQUIREMENTS,
    entry_points={'console_scripts': ['postqe=postqe.cli:main']},
    cmdclass={
        'install_lib': InstallLibCommand,
        'clean': CleanCommand,
    },
    author="Mauro Palumbo, Pietro Delugas, Davide Brunato",
    author_email="pdelugas@sissa.it, brunato@sissa.it",
    license='LGPL-2.1',
    long_description="Post processing tools for Quantum Espresso",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
