from setuptools import setup, find_packages
from distutils.core import Extension
from Cython.Build import cythonize
from majiq.src.constants import VERSION
import numpy
try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')


extensions = [Extension('majiq.src.normalize', ['majiq/src/normalize.pyx'])]
extensions += [Extension('majiq.src.plotting', ['majiq/src/plotting.pyx'])]
extensions += [Extension('majiq.src.polyfitnb', ['majiq/src/polyfitnb.pyx'])]
extensions += [Extension('majiq.src.sample', ['majiq/src/sample.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.psi', ['majiq/src/psi.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.analize', ['majiq/src/analize.pyx'], include_dirs=[numpy.get_include()])]
extensions += [Extension('majiq.src.io_base', ['majiq/src/io_base.pyx'], include_dirs=pysam.get_include())]

include_dirs=pysam.get_include()

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
