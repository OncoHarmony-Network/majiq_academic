import os
import sys
import numpy as np

from setuptools import Extension, setup


# for htslib support
HTSLIB_LIBRARY = ["hts", "z"]
HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", "/usr/local/include")]
HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", "/usr/local/lib")]

# numpy include directories
NPY_INC_DIRS = [np.get_include()]
# currently using deprecated API (since numpy v1.7)
# silence deprecation warnings for now
SILENCE_DEPRECATION = [("NPY_NO_DEPRECATED_API", 0)]

# majiq include/library directories
MAJIQ_INC_DIRS = ["rna_majiq/src/internals", "rna_majiq/c"]
MAJIQ_INC_STATS_DIR = ["rna_majiq/src/internals/stats"]
MAJIQ_LIB_DIRS = ["rna_majiq/src/internals"]

# openmp support flags
compile_args = ["-fopenmp"]
linker_args = ["-lgomp"]

# uncomment to compile with symbols for debugging with gdb/lldb
# compile_args += ["-O0", "-g"]  # for debugging

# C++11 flags
if sys.platform == "darwin":
    pass
else:
    compile_args.append("-std=c++11")


# define Cython extension modules
ext_modules = [
    Extension(
        "rna_majiq.src.build",
        [
            "rna_majiq/src/build.pyx",
            "rna_majiq/src/internals/io_bam.cpp",
            "rna_majiq/src/internals/grimoire.cpp",
        ],
        include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS + HTSLIB_INC_DIRS,
        library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
        libraries=HTSLIB_LIBRARY,
        runtime_library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
        define_macros=SILENCE_DEPRECATION,
        language="c++",
    ),
    Extension(
        "rna_majiq.src.calc_psi",
        ["rna_majiq/src/calc_psi.pyx", "rna_majiq/src/internals/psi.cpp"],
        include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
        library_dirs=MAJIQ_LIB_DIRS,
        runtime_library_dirs=MAJIQ_LIB_DIRS,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
        define_macros=SILENCE_DEPRECATION,
        language="c++",
    ),
    Extension(
        "rna_majiq.src.deltapsi",
        ["rna_majiq/src/deltapsi.pyx", "rna_majiq/src/internals/psi.cpp"],
        include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
        library_dirs=MAJIQ_LIB_DIRS,
        runtime_library_dirs=MAJIQ_LIB_DIRS,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
        define_macros=SILENCE_DEPRECATION,
        language="c++",
    ),
    Extension(
        "rna_majiq.src.indpnt",
        ["rna_majiq/src/indpnt.pyx", "rna_majiq/src/internals/psi.cpp"],
        include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS + MAJIQ_INC_STATS_DIR,
        library_dirs=MAJIQ_LIB_DIRS,
        runtime_library_dirs=MAJIQ_LIB_DIRS,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
        define_macros=SILENCE_DEPRECATION,
        language="c++",
    ),
    Extension(
        "rna_majiq.src.io",
        ["rna_majiq/src/io.pyx"],
        include_dirs=NPY_INC_DIRS + MAJIQ_INC_DIRS,
        define_macros=SILENCE_DEPRECATION,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
        language="c++",
    ),
    Extension(
        "rna_majiq.c.splice_graph_sql",
        ["rna_majiq/c/splice_graph_sql.pyx", "rna_majiq/c/sqlite3.c"],
        language="c++",
        include_dirs=NPY_INC_DIRS,
        extra_compile_args=compile_args,
        extra_link_args=linker_args,
    ),
]


for x in ext_modules:
    x.cython_directives = {"language_level": "3"}  # Python3, not 2


setup(
    ext_modules=ext_modules,
    use_scm_version={
        "root": "..",
        "relative_to": __file__,
        "fallback_version": "2.5.5+scmfallback",
    },
)
