#
# Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the LGPL-2.1 license. See the
# file 'LICENSE' in the root directory of the present distribution, or
# https://opensource.org/licenses/LGPL-2.1
#

import os
import shutil
import re
import platform
from pathlib import Path

from setuptools import setup, Distribution
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

from distutils import log
from distutils.command.clean import clean  # type: ignore[attr-defined]
from distutils.file_util import copy_file


VERSION_NUMBER_PATTERN = re.compile(r"version_number\s*=\s*(\'[^\']*\'|\"[^\"]*\")")

# QE_SOURCE_URL = "https://github.com/QEF/q-e/archive/refs/tags/qe-7.2.zip"
# QE_SOURCE_MD5SUM = "787c4aad3b203f7dd61d03a849a1c4e9"

with Path(__file__).parent.joinpath("requirements.txt").open() as fp:
    REQUIREMENTS = fp.read()


def find_extension_module(pattern):
    """
    Returns the absolute pathname of an extension module built for the running platform.

    :param: a glob pattern for finding the candidate extension modules.
    :return: A pathname string or `None` if no suitable module is found.
    """
    python_version = "".join(platform.python_version_tuple()[:2])

    for filename in map(str, Path(__file__).parent.glob(pattern)):
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


class BuildExtCommand(build_ext):
    """
    BuildExtCommand is a custom extension of the build_ext command for setuptools.

    This class is responsible for building the Fortran modules and Quantum Espresso library,
    as well as the pyqe Python module, during the package build process.

    The following steps are performed during the build process:

    - Prepare the build directory.
    - Copy the Fortran source from the package directory to the build directory.
    - Build f90utils Fortran module.
    - Build Quantum Espresso library.
    - Build pyqe Python module.

    Attributes: build_temp (str): The temporary build directory.
                inplace (bool): If True, the build is done in-place.
                verbose (bool): If True, additional information is printed during build process.
                dry_run (bool): If True, the build process will not actually be performed.

    Methods: run(): Perform the build process.
             _build_f90utils_module(build_dir: Path): Build the f90utils Fortran module.
             _build_quantum_espresso(qe_build_dir: Path): Build Quantum Espresso library.
             _build_pyqe_module(build_dir: Path): Build pyqe Python module.
    """

    def run(self):
        build_ext.run(self)

        build_dir = Path(self.build_temp).absolute().joinpath("postqe-fortran")
        qe_build_dir = build_dir.joinpath("qe")
        if not qe_build_dir.exists():
            qe_build_dir.mkdir(parents=True)

        build_py = self.get_finalized_command("build_py")
        package_dir = Path(build_py.get_package_dir("postqe"))
        package_fortran_dir = str(package_dir.joinpath("fortran"))

        for item in os.listdir(package_fortran_dir):
            # Copy everything its inside the package postqe/fortran folder into
            # the postqe-fortran folder inside the build dicrectory
            source_item = os.path.join(package_fortran_dir, item)
            destination_item = os.path.join(build_dir, item)

            if os.path.isdir(source_item):
                if os.path.exists(destination_item):
                    shutil.rmtree(destination_item)
                shutil.copytree(source_item, destination_item)
            else:
                shutil.copy2(source_item, destination_item)

        self._build_f90utils_module(build_dir)
        self._build_quantum_espresso(qe_build_dir)
        self._build_pyqe_module(build_dir)

        if self.inplace:
            copy_file(
                src=str(build_dir.joinpath("pyqe.py")),
                dst=str(package_dir),
                verbose=self.verbose,
                dry_run=self.dry_run,
            )
            for filename in build_dir.glob("_pyqe*.so"):
                copy_file(
                    src=str(filename),
                    dst=str(package_dir),
                    verbose=self.verbose,
                    dry_run=self.dry_run,
                )
            for filename in build_dir.glob("f90utils*.so"):
                copy_file(
                    src=str(filename),
                    dst=str(package_dir),
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

        # Check QE installation
        for version_file in qe_topdir.glob("include/*version.h"):
            with version_file.open() as fp:
                version_number = VERSION_NUMBER_PATTERN.search(fp.read())

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
        qe_build_dir = build_dir.joinpath("qe")
        # Read the CMakeCache.txt file to find:
        # LAPACK LIBS path
        # BLAS LIBS path
        # FFTW3 LIBS path
        with open(f"{qe_build_dir}/CMakeCache.txt", "r") as file:
            cmakecache_lines = file.readlines()

        pattern = r"(/[^;\s\]]+)"

        fftw3_first_path = None
        lapack_path = None
        blas_path = None

        for line in cmakecache_lines:
            if (
                line.startswith("FIND_PACKAGE_MESSAGE_DETAILS_FFTW3")
                and fftw3_first_path is None
            ):
                match = re.search(pattern, line)
                if match:
                    fftw3_first_path = match.group(1)

            if (
                line.startswith("FIND_PACKAGE_MESSAGE_DETAILS_LAPACK")
                and lapack_path is None
                and blas_path is None
            ):
                matches = re.findall(pattern, line)
                if len(matches) == 2:
                    lapack_path, blas_path = matches

            # Break the loop if both paths are found
            if fftw3_first_path and lapack_path and blas_path:
                break

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
            os.environ["BLAS_LIBS"] = blas_path
        else:
            print("No BLAS path found")

        # Run the make target to build the pyqe module
        os.system(
            f"make LAPACK_LIBS={os.environ['LAPACK_LIBS']} "
            f"BLAS_LIBS={os.environ['BLAS_LIBS']} "
            f"FFTW3_LIBS={os.environ['FFTW3_LIBS']} "
            f"QE_BUILD_DIR={qe_build_dir} "
            f"BUILD_DIR={build_dir} "
            f"-f pyqe_Makefile "
            f"-C {build_dir} "
            f"pyqe"
        )
        # Modify python wrapper module, fixing import in postqe package.
        with build_dir.joinpath("pyqe.py").open() as fp:
            python_wrapper_lines = fp.readlines()

        with build_dir.joinpath("pyqe.py").open(mode="w") as fp:
            python_wrapper_lines[0] = "# Altered wrapper for postqe\n"
            python_wrapper_lines[1] = "from . import _pyqe\n"
            fp.writelines(python_wrapper_lines)


class InstallCommand(install):
    distribution: Distribution  # for avoid static check warning

    def run(self):
        if find_extension_module("postqe/_pyqe*.so") is None or find_extension_module(
            "postqe/f90utils*.so"
        ):
            print("A required extension module not found, invoke build_ext ...")
            cmd_obj = self.distribution.get_command_obj("build_ext")
            cmd_obj.inplace = True
            self.run_command("build_ext")
        install.run(self)


class CleanCommand(clean):
    def run(self):
        clean.run(self)
        if self.all:
            log.info("Remove built extension modules ...")
            source_dir = Path(__file__).parent.joinpath("postqe")

            wrapper_module = source_dir.joinpath("pyqe.py")
            if wrapper_module.is_file():
                log.info("Remove {}".format(str(wrapper_module)))
                wrapper_module.unlink()

            for ext_module in source_dir.glob("*.so"):
                log.info(f"Remove {ext_module}")
                os.unlink(ext_module)


setup(
    name="postqe",
    version="1.0.0",
    packages=["postqe"],
    package_data={"postqe": ["_pyqe.*.so", "f90utils*.so"]},
    install_requires=REQUIREMENTS,
    entry_points={"console_scripts": ["postqe=postqe.cli:main"]},
    cmdclass={
        "build_ext": BuildExtCommand,
        "install": InstallCommand,
        "clean": CleanCommand,
    },
    author="Mauro Palumbo, Pietro Delugas, Davide Brunato",
    author_email="pdelugas@sissa.it, brunato@sissa.it",
    license="LGPL-2.1",
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
