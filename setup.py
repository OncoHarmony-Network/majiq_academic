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

compile_args = ['-g', '-std=c++11']
linker_args = ['-std=c++11']
if sys.platform == 'darwin':
    compile_args.append('-stdlib=libc++')
    linker_args = ['-L/usr/local/opt/llvm/lib']

pysam_library_path = [os.path.abspath(os.path.join(os.path.dirname(pysam.__file__)))]
extensions = [Extension('majiq.src.internals.seq_parse',
                        ['majiq/src/internals/seq_parse.pyx', 'majiq/src/internals/io_bam.cpp'],
                        include_dirs=['majiq/src/internals'] + pysam.get_include(),
                        library_dirs=pysam_library_path,
                        libraries=['htslib'],
                        runtime_library_dirs=pysam_library_path,
                        extra_compile_args=compile_args,  extra_link_args=linker_args,
                        language='c++')]


extensions += [Extension('majiq.src.normalize', ['majiq/src/normalize.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.grimoire.junction', ['majiq/grimoire/junction.pyx'])]
extensions += [Extension('majiq.grimoire.lsv', ['majiq/grimoire/lsv.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.grimoire.exon', ['majiq/grimoire/exon.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.plotting', ['majiq/src/plotting.pyx'])]
extensions += [Extension('majiq.src.polyfitnb', ['majiq/src/polyfitnb.pyx'], include_dirs=[numpy.get_include()])]
inc_dirs = [numpy.get_include()]
inc_dirs.extend(pysam.get_include())
extensions += [Extension('majiq.src.io_bam', ['majiq/src/io_bam.pyx'], include_dirs=inc_dirs)]
extensions += [Extension('majiq.src.io', ['majiq/src/io.pyx'], include_dirs=inc_dirs)]
extensions += [Extension('majiq.src.sample', ['majiq/src/sample.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.adjustdelta', ['majiq/src/adjustdelta.pyx'], include_dirs=[numpy.get_include()])]


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
