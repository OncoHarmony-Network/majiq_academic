from setuptools import setup
from rna_voila.constants import VERSION, store_git_version
import os

script_path = os.path.dirname(__file__)
if script_path:
    os.chdir(script_path)

try:
    from git import Repo

    repo = Repo('./')
    assert not repo.bare
    sha = repo.head.object.hexsha
    short_sha = repo.git.rev_parse(sha, short=7)
    store_git_version(short_sha)
except Exception as e:
    print(f"Problem was encounted during git hash extraction. Hash value disabled: {e}")


def requirements():
    with open('requirements.txt', 'r') as f:
        reql = [l.strip() for l in f.readlines()]
        print(reql)
        return reql


def files_recursive(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


setup(
    name='rna_voila',
    version=VERSION,
    description="VOILA -- visualization tool for MAJIQ",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    license='No License',
    packages=[
        'rna_voila',
        'rna_voila.api',
        'rna_voila.view',
        'rna_voila.utils',
        'rna_voila.view',
    ],
    install_requires=requirements(),
    include_package_data=True,
    package_data={
        'rna_voila':
        ['api/model.sql']
        + files_recursive('rna_voila/view/templates')
        + files_recursive('rna_voila/view/static')
    },
    zip_safe=False,
    entry_points={'console_scripts': ['voila = rna_voila.run_voila:main']},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Bioinformaticians',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.8']
)
