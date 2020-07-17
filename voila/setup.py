from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools import Extension
from rna_voila.constants import VERSION, store_git_version
import sys
import os

script_path = os.path.dirname(__file__)
os.chdir(script_path)

try:
    from git import Repo

    repo = Repo('./')
    assert not repo.bare
    sha = repo.head.object.hexsha
    short_sha = repo.git.rev_parse(sha, short=7)
    store_git_version(short_sha)

except Exception as e:
    print ("Problem was encounted during it hash extraction. Hash value disabled: %s" % e)


def requirements():
    with open('requirements.txt', 'r') as f:
        reql = [l.strip() for l in f.readlines()]
        print(reql)
        return reql


def setup_requirements(search=["Cython", "numpy"]):
    return [x for x in requirements() if any(y in x for y in search)] + ['wheel']


class InstallCommand(install):

    user_options = install.user_options + [
        ("num-threads=", "j", "Max number of threads to use for compiling extensions")
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.num_threads = 0

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        print(self.__dict__)

        self.distribution.packages = ['rna_voila', 'rna_voila.api', 'rna_voila.view', 'rna_voila.utils',
                                      'rna_voila.view']

        install.run(self)

def files_recursive(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

setup(
    name = 'rna_voila',
    version = VERSION,
    description="VOILA -- visualization tool for MAJIQ",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    license='No License',
    packages=find_packages(),
    setup_requires=setup_requirements(),
    install_requires=requirements(),
    include_package_data=True,
    package_data={'rna_voila': ['api/model.sql'] +
                           files_recursive('rna_voila/view/templates') +
                           files_recursive('rna_voila/view/static')},
    zip_safe=False,
    entry_points={'console_scripts': ['voila = rna_voila.run_voila:main']},
    cmdclass={'install': InstallCommand},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Bioinformaticians',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.8']
)

