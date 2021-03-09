import os
from setuptools import setup

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension
import pybind11.setup_helpers as sh


__version__ = "0.0.1"


# disable pybind11 helpers for OSX -- they are not compiler independent
sh.MACOS = False

# for htslib support
HTSLIB_LIBRARY = ["hts", "z"]
HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", "/usr/local/include")]
HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", "/usr/local/lib")]

ext_modules = [
    Pybind11Extension(
        "new_majiq.internals",
        sorted([
            "src/internals/PassedJunctions.cpp",
            "src/internals/SpliceGraph.cpp",
            "src/internals/ContigIntrons.cpp",
            "src/internals/SJJunctionsPositions.cpp",
            "src/internals/GFF3.cpp",
            "src/internals/TranscriptModels.cpp",
            "src/pySpliceGraph.cpp",
            "src/new_majiq.cpp",
        ]),
        define_macros=[('VERSION_INFO', __version__), ("DEBUG", )],
        include_dirs=["src/", *HTSLIB_INC_DIRS],
        library_dirs=[*HTSLIB_LIB_DIRS],
        libraries=[*HTSLIB_LIBRARY],
        cxx_std=17,
        # extra_compile_args=["-O0", "-g"],
    ),
]


setup(
    version=__version__,
    packages=["new_majiq"],
    package_dir={"": "src"},
    ext_modules=ext_modules,
)
