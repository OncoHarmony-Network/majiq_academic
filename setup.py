from setuptools import setup, find_packages
from distutils.core import Extension
from Cython.Build import cythonize
from majiq.src.constants import VERSION
import numpy
import sys
import os

print("COMPILER ", os.environ.get("CXX", None))

extensions = []
HTSLIB_LIBRARY = ['hts', 'z']
HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", '/usr/local/lib')]
HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", '/usr/local/include')]


compile_args = ['-fopenmp']
scythe_compiler_args = ['-DSCYTHE_COMPILE_DIRECT', '-DSCYTHE_PTHREAD']
linker_args = ['-lgomp']

if sys.platform == 'darwin':
    # pass
    compile_args.append('-stdlib=libc++')
    linker_args = ['-L/usr/local/opt/llvm/lib'] + linker_args
else:
    compile_args.append('-std=c++11')


MAJIQ_INC_DIRS = ['majiq/src/internals']
VOILA_INC_DIRS = ['voila/c']
MAJIQ_LIB_DIRS = ['majiq/src/internals']
NPY_INC_DIRS   = [numpy.get_include()]

extensions += [Extension('majiq.src.polyfitnb', ['majiq/src/polyfitnb.pyx'], language='c++', include_dirs=NPY_INC_DIRS)]
extensions += [Extension('majiq.src.internals.seq_parse', ['majiq/src/internals/seq_parse.pyx',
                                                           'majiq/src/internals/io_bam.cpp',
                                                           'majiq/src/internals/grimoire.cpp'],
                         include_dirs=MAJIQ_INC_DIRS + VOILA_INC_DIRS + NPY_INC_DIRS + HTSLIB_INC_DIRS,
                         library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                         libraries=HTSLIB_LIBRARY,
                         runtime_library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                         extra_compile_args=compile_args+scythe_compiler_args,  extra_link_args=linker_args,
                         language='c++')]

extensions += [Extension('majiq.src.calc_psi', ['majiq/src/calc_psi.pyx',
                                                'majiq/src/internals/psi.cpp'],
                         include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                         library_dirs=MAJIQ_LIB_DIRS,
                         runtime_library_dirs=MAJIQ_LIB_DIRS,
                         extra_compile_args=compile_args+scythe_compiler_args,  extra_link_args=linker_args,
                         language='c++')]

extensions += [Extension('majiq.src.deltapsi', ['majiq/src/deltapsi.pyx',
                                                'majiq/src/internals/psi.cpp'],
                         include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                         library_dirs=MAJIQ_LIB_DIRS,
                         runtime_library_dirs=MAJIQ_LIB_DIRS,
                         extra_compile_args=compile_args+scythe_compiler_args,  extra_link_args=linker_args,
                         language='c++')]

extensions += [Extension('majiq.src.io', ['majiq/src/io.pyx'], language='c++',
                         include_dirs=NPY_INC_DIRS + MAJIQ_INC_DIRS,
                         extra_compile_args=compile_args,  extra_link_args=linker_args,)]

# extensions += [Extension('majiq.src.normalize', ['majiq/src/normalize.pyx'],     include_dirs=NPY_INC_DIRS)]
extensions += [Extension('majiq.src.plotting', ['majiq/src/plotting.pyx'])]
extensions += [Extension('majiq.src.psi', ['majiq/src/psi.pyx'],                 include_dirs=NPY_INC_DIRS)]
extensions += [Extension('majiq.src.sample', ['majiq/src/sample.pyx'],           include_dirs=NPY_INC_DIRS)]
extensions += [Extension('majiq.src.adjustdelta', ['majiq/src/adjustdelta.pyx'], include_dirs=NPY_INC_DIRS)]
extensions += [Extension('voila.c.splice_graph_sql', ['voila/c/splice_graph_sql.pyx', 'voila/c/sqlite3.c'],
                         language='c++', include_dirs=NPY_INC_DIRS, extra_compile_args=compile_args,
                         extra_link_args=linker_args)]

setup(
    name="majiq",
    packages=find_packages(),
    version=VERSION,
    description="MAJIQ and VOILA",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    license='LICENSE.txt',
    include_package_data=True,
    entry_points={'console_scripts': ['majiq = majiq.run_majiq:main', 'voila = voila.run_voila:main']},
    zip_safe=False,
    ext_modules=cythonize(extensions),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Bioinformaticians',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.4']
)
