import numpy
from Cython.Build import cythonize
from setuptools import setup, find_packages, Extension

from majiq.src.constants import VERSION

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
    entry_points={'console_scripts': ['majiq = majiq.run_majiq:main', 'voila = voila.run_voila:main']},
    zip_safe=False,
    ext_modules=cythonize(extensions),
    install_requires=[
        'appdirs==1.4.3',
        'click==6.7',
        'colorama==0.3.9',
        'cycler==0.10.0',
        'Cython==0.25.2',
        'decorator==4.0.11',
        'h5py==2.7.0',
        'Jinja2==2.9.6',
        'MarkupSafe==1.0',
        'matplotlib==2.0.2',
        'networkx==1.11',
        'numpy==1.13.1',
        'packaging==16.8',
        'pyparsing==2.2.0',
        'pysam==0.12',
        'python-dateutil==2.6.0',
        'pytz==2017.2',
        'scipy==0.19.0',
        'six==1.10.0',
        'pandas==0.20.2',
        'quicksect==0.0.2',
        'SQLAlchemy==1.2.1',
        'psutil==5.4.3',
    ],
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
