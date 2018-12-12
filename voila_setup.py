#!/usr/bin/env python
"""
This setup file can be used to install only voila.  This can be used with the following command:

    python setup_voila.py install

"""
from setuptools import setup

from voila.constants import VERSION

setup(
    name='voila',
    version=VERSION,
    packages=['voila', 'voila.api', 'voila.view', 'voila.utils', 'voila.view'],
    install_requires=['Flask-WTF', 'Flask', 'h5py', 'scipy', 'waitress'],
    entry_points={'console_scripts': ['voila = voila.run_voila:main']},
    package_data={
        'voila': [
            'view/templates/*',
            'view/static/css/*',
            'view/static/js/*',
            'view/static/img/*'
        ]
    }
)
