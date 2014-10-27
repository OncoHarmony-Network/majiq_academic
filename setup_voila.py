#!/usr/bin/env python

from cx_Freeze import setup, Executable

packages =['jinja2', 'numpy']
include_files = [('voila/templates', 'templates')]
excludes = ['matplotlib']
setup(name='Voila',
      version='0.8.0',
      description='Visualization of Local Splice Variants',
      author='Alejandro Barrera',
      author_email='abarrera@biociphers.org',
      url='http://www.biociphers.org/',
      options = { 'build_exe':
                      {'packages': packages,
                       'include_files': include_files,
                       'excludes': excludes
                      },
                  # 'iconfile': ['/Users/abarrera/majiq/voila/templates/static/img/icons/biociphers.ico']
                },
      executables= [Executable(script='voila/run_voila.py', targetName='voila')]
     )

# python setup_voila.py build_exe