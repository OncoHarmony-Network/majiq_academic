import pip

pip.main(['install', 'Cython==0.25.2'])
pip.main(['install', '-r', 'requirements.txt', '-U'])

from setuptools import setup, find_packages
from distutils.extension import Extension
from majiq.src.constants import VERSION
from Cython.Distutils import build_ext
import numpy

try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')

extensions = [Extension('majiq.src.normalize', ['majiq/src/normalize.pyx'], include_dirs=[numpy.get_include()])]
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
extensions += [Extension('majiq.src.psi', ['majiq/src/psi.pyx'], include_dirs=[numpy.get_include()])]

include_dirs = pysam.get_include()

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
    cmdclass={'build_ext': build_ext},
    ext_modules=extensions,
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
