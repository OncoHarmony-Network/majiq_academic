#!/usr/bin/env python
from setuptools import setup
from voila.constants import VERSION

setup(
    name='voila',
    version=VERSION,
    packages=['voila', 'voila.api', 'voila.view', 'voila.utils', 'voila.view'],
    install_requires=['h5py', 'flask', 'sqlalchemy', 'waitress', 'scipy'],
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
