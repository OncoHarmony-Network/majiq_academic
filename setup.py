from setuptools import setup, find_packages
from distutils.core import Extension
from Cython.Build import cythonize
from majiq.src.constants import VERSION
import numpy
import sys
import os
try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')

compile_args = ['-fopenmp', '-std=c++11']
linker_args = ['-lgomp', '-std=c++11']

if sys.platform == 'darwin':
    compile_args.append('-stdlib=libc++')
    linker_args = ['-L/usr/local/opt/llvm/lib'] + linker_args
# #
# os.environ['CC'] = 'g++-8'
# os.environ['CXX'] = 'g++-8'

include_librs = ['majiq/src/internals', numpy.get_include()] + pysam.get_include()
pysam_library_path = [os.path.abspath(os.path.join(os.path.dirname(pysam.__file__)))]
extensions = [Extension('majiq.src.internals.seq_parse',
                        ['majiq/src/internals/seq_parse.pyx', 'majiq/src/internals/io_bam.cpp',
                         'majiq/src/internals/grimoire.cpp'],
                        include_dirs=include_librs,
                        library_dirs=pysam_library_path,
                        libraries=['htslib'],
                        runtime_library_dirs=pysam_library_path,
                        extra_compile_args=compile_args,  extra_link_args=linker_args,
                        language='c++', gdb_debug=True)]

extensions += [Extension('majiq.src.io', ['majiq/src/io.pyx'], language='c++', include_dirs=include_librs,
                         extra_compile_args=compile_args,  extra_link_args=linker_args,)]

extensions += [Extension('majiq.src.normalize', ['majiq/src/normalize.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.grimoire.junction', ['majiq/grimoire/junction.pyx'])]
extensions += [Extension('majiq.grimoire.lsv', ['majiq/grimoire/lsv.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.grimoire.exon', ['majiq/grimoire/exon.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.plotting', ['majiq/src/plotting.pyx'])]
extensions += [Extension('majiq.src.polyfitnb', ['majiq/src/polyfitnb.pyx'], include_dirs=[numpy.get_include()])]
inc_dirs = [numpy.get_include()]
inc_dirs.extend(pysam.get_include())
extensions += [Extension('majiq.src.psi', ['majiq/src/psi.pyx'], include_dirs=inc_dirs)]
extensions += [Extension('majiq.src.sample', ['majiq/src/sample.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.adjustdelta', ['majiq/src/adjustdelta.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('voila.c.splice_graph_sql', ['voila/c/splice_graph_sql.pyx'])]

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
