#!/usr/bin/env python

from cx_Freeze import setup, Executable

packages = ['scipy', 'pysam', 'numpy', 'matplotlib']
include_files = []
excludes = []
setup(name='Majiq',
      version='0.8.0',
      description='MAJIQ',
      author='Jordi Vaquero',
      author_email='jordi@biociphers.org',
      url='http://www.biociphers.org/',
      options={'build_exe': {'packages': packages,
                             'include_files': include_files,
                             'excludes': excludes}
               },
      executables=[Executable(script='majiq')])

# python setup_majiq.py build_exe