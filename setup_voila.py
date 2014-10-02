#!/usr/bin/env python

from cx_Freeze import setup, Executable

packages =['jinja2']
include_files = [('voila/templates', 'templates')]
excludes = ['matplotlib', 'scipy']
setup(name='Voila',
      version='0.1.0',
      description='Visualization of Local Splice Variants',
      author='Alejandro Barrera',
      author_email='abarrera@biociphers.org',
      url='http://www.biociphers.org/',
      options = { 'build_exe':
                      {'packages': packages,
                       'include_files': include_files,
                       'excludes': excludes
                      },
                  'iconfile': '/Users/abarrera/majiq/voila/templates/static/img/icons/biociphers.ico'
                },
      executables= [Executable('voila/run_voila.py')]
     )