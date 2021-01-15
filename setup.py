from setuptools import setup

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension
import pybind11.setup_helpers as sh


__version__ = "0.0.1"


# disable pybind11 helpers for OSX -- they are not compiler independent
sh.MACOS = False


ext_modules = [
    Pybind11Extension(
        "new_majiq",
        sorted([
            "src/pySpliceGraph.cpp",
            "src/new_majiq.cpp",
        ]),
        define_macros=[('VERSION_INFO', __version__), ("DEBUG", )],
        include_dirs=["src/"],
        cxx_std=17,
    ),
]


setup(
    name="new_majiq",
    version=__version__,
    author="Joseph K Aicher",
    description="Python bindings into new_majiq c++ code",
    packages=[],
    package_dir={"": "src"},
    ext_modules=ext_modules,
    zip_safe=False,
)
