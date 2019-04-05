from setuptools import setup, find_packages
from setuptools.command.install import install
from distutils.core import Extension
from Cython.Build import cythonize
from majiq.src.constants import VERSION, store_git_version
import numpy
import sys
import os

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
    with open('requirements.txt') as f:
        return [l.strip() for l in f.readlines()]

class InstallCommand(install):

    user_options = install.user_options + [('voila-only', 'v', None),]

    def initialize_options(self):
        install.initialize_options(self)
        self.voila_only = 0

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        print(self.voila_only)
        if(self.voila_only):
            # self.name = 'voila',
            # self.version = VERSION,
            self.distribution.packages = ['voila', 'voila.api', 'voila.view', 'voila.utils', 'voila.view']
            self.distribution.package_data = {'voila': ['view/templates/*', 'view/static/css/*',
                                           'view/static/js/*', 'view/static/img/*']}
        else:
            extensions = []
            HTSLIB_LIBRARY = ['hts', 'z']
            HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", '/usr/local/lib')]
            HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", '/usr/local/include')]

            compile_args = ['-fopenmp']
            scythe_compiler_args = ['-DSCYTHE_COMPILE_DIRECT', '-DSCYTHE_PTHREAD']
            linker_args = ['-lgomp']

            if sys.platform == 'darwin':
                # os.environ['CLANG_DEFAULT_CXX_STDLIB'] = 'libc++'
                pass
                # compile_args.append('-stdlib=libc++')
            #     linker_args = ['-L/usr/local/opt/llvm/lib'] + linker_args
            else:
                compile_args.append('-std=c++11')

            MAJIQ_INC_DIRS = ['majiq/src/internals']
            MAJIQ_INC_STATS_DIR = ['majiq/src/internals/stats']
            VOILA_INC_DIRS = ['voila/c']
            MAJIQ_LIB_DIRS = ['majiq/src/internals']

            NPY_INC_DIRS = [numpy.get_include()]

            extensions += [Extension('majiq.src.polyfitnb', ['majiq/src/polyfitnb.pyx'], language='c++',
                                     include_dirs=NPY_INC_DIRS)]
            extensions += [Extension('majiq.src.build', ['majiq/src/build.pyx',
                                                         'majiq/src/internals/io_bam.cpp',
                                                         'majiq/src/internals/grimoire.cpp'],
                                     include_dirs=MAJIQ_INC_DIRS + VOILA_INC_DIRS + NPY_INC_DIRS + HTSLIB_INC_DIRS,
                                     library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                                     libraries=HTSLIB_LIBRARY,
                                     runtime_library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                                     extra_compile_args=compile_args + scythe_compiler_args,
                                     extra_link_args=linker_args,
                                     language='c++')]

            extensions += [Extension('majiq.src.calc_psi', ['majiq/src/calc_psi.pyx', 'majiq/src/internals/psi.cpp'],
                                     include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                                     library_dirs=MAJIQ_LIB_DIRS,
                                     runtime_library_dirs=MAJIQ_LIB_DIRS,
                                     extra_compile_args=compile_args + scythe_compiler_args,
                                     extra_link_args=linker_args,
                                     language='c++')]

            extensions += [Extension('majiq.src.deltapsi', ['majiq/src/deltapsi.pyx', 'majiq/src/internals/psi.cpp'],
                                     include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                                     library_dirs=MAJIQ_LIB_DIRS,
                                     runtime_library_dirs=MAJIQ_LIB_DIRS,
                                     extra_compile_args=compile_args + scythe_compiler_args,
                                     extra_link_args=linker_args,
                                     language='c++')]

            extensions += [Extension('majiq.src.indpnt', ['majiq/src/indpnt.pyx', 'majiq/src/internals/psi.cpp'],
                                     include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS + MAJIQ_INC_STATS_DIR,
                                     library_dirs=MAJIQ_LIB_DIRS,
                                     runtime_library_dirs=MAJIQ_LIB_DIRS,
                                     extra_compile_args=compile_args + scythe_compiler_args,
                                     extra_link_args=linker_args,
                                     language='c++')]

            extensions += [Extension('majiq.src.io', ['majiq/src/io.pyx'], include_dirs=NPY_INC_DIRS + MAJIQ_INC_DIRS,
                                     extra_compile_args=compile_args, extra_link_args=linker_args, language='c++')]

            extensions += [Extension('voila.c.splice_graph_sql', ['voila/c/splice_graph_sql.pyx', 'voila/c/sqlite3.c'],
                                     language='c++', include_dirs=NPY_INC_DIRS, extra_compile_args=compile_args,
                                     extra_link_args=linker_args)]

            self.distribution.packages = find_packages()
            self.distribution.include_package_data = True
            self.distribution.entry_points['console_scripts'].append('majiq = majiq.run_majiq:main')
            self.distribution.zip_safe = False
            self.distribution.ext_modules = cythonize(extensions, language_level=3)


        install.run(self)

setup(
    name = 'majiq',
    version = VERSION,
    description="MAJIQ and VOILA",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    license='LICENSE.txt',
    entry_points = {'console_scripts': ['voila = voila.run_voila:main']},
    install_requires = requirements(),
    cmdclass={'install': InstallCommand,},
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

