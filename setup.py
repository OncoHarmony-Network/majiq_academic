import os

# get include path for python distribution
from sysconfig import get_paths

import numpy as np  # and for numpy
import pybind11.setup_helpers as sh

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension
from setuptools import Extension, setup

cpp_version = "0.0.1"


# disable pybind11 helpers for OSX -- they are not compiler independent
sh.MACOS = False

# for htslib support
HTSLIB_LIBRARY = ["hts", "z"]
HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", "/usr/local/include")]
HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", "/usr/local/lib")]

ext_modules = [
    Pybind11Extension(
        "new_majiq.internals",
        sorted(
            [
                "src/internals/PassedJunctions.cpp",
                "src/internals/SpliceGraph.cpp",
                "src/internals/SJIntrons.cpp",
                "src/internals/SJBinsReads.cpp",
                "src/internals/GFF3.cpp",
                "src/internals/TranscriptModels.cpp",
                "src/pySpliceGraph.cpp",
                "src/new_majiq.cpp",
            ]
        ),
        define_macros=[("VERSION_INFO", cpp_version), ("DEBUG",)],
        include_dirs=["src/", *HTSLIB_INC_DIRS],
        library_dirs=[*HTSLIB_LIB_DIRS],
        runtime_library_dirs=[*HTSLIB_LIB_DIRS],
        libraries=[*HTSLIB_LIBRARY],
        cxx_std=17,
        # extra_compile_args=["-O0", "-g"],
    ),
    Extension(
        "new_majiq.beta_mixture",
        sources=["src/gufuncs/beta_mixture/beta_mixture_module.cpp"],
        language="c++",
        include_dirs=[
            "src/",
            get_paths()["include"],
            np.get_include(),
        ],
        # extra_compile_args=["-g0", "-std=c++11"],
        extra_compile_args=["-g0", "-std=c++11"],
        extra_link_args=["-std=c++11"],
    ),
    Extension(
        "new_majiq._offsets",
        sources=["src/gufuncs/offsets/offsets_module.cpp"],
        language="c++",
        include_dirs=[
            "src/",
            get_paths()["include"],
            np.get_include(),
        ],
        extra_compile_args=["-g0", "-std=c++17"],
        extra_link_args=["-std=c++17"],
    ),
    Extension(
        "new_majiq._stats",
        sources=["src/gufuncs/stats/stats_module.cpp"],
        language="c++",
        include_dirs=[
            "src/",
            get_paths()["include"],
            np.get_include(),
        ],
        extra_compile_args=["-g0", "-std=c++17"],
        extra_link_args=["-std=c++17"],
    ),
]


setup(
    ext_modules=ext_modules,
)
