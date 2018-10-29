#!/usr/bin/env python
from setuptools import setup
from voila.constants import VERSION

setup(
    name='voila',
    version=VERSION,
    packages=['voila', 'voila.api', 'voila.flask_proj', 'voila.utils', 'voila.view'],
    install_requires=['h5py', 'flask', 'sqlalchemy', 'waitress', 'scipy'],
    entry_points={'console_scripts': ['voila = voila.run_voila:main']},
    package_data={
        'voila': [
            'flask_proj/templates/*',
            'flask_proj/static/css/*',
            'flask_proj/static/js/*',
            'flask_proj/static/img/*'
        ]
    }
)
