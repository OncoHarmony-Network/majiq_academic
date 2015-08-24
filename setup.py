from setuptools import setup, find_packages

setup(
    name="majiq",
    packages=find_packages('majiq'),
    version="0.8.0",
    description="MAJIQ",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    include_package_data=True,
    entry_points={'console_scripts': ['majiq = majiq.parser:main', 'voila = majiq.voila.run_voila:main']},
    zip_safe=False,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Bioinformaticians',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
    ]
)
