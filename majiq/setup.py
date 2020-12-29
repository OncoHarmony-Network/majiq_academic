from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools import Extension
from rna_majiq.src.constants import VERSION, store_git_version
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
    with open('requirements.txt') as f:
        reql = [l.strip() for l in f.readlines()]
        # print(reql)
        return reql


def setup_requirements(search=["Cython", "numpy", "voila"]):
    return [x for x in requirements() if any(y in x for y in search)] + ['wheel']


class InstallCommand(install):

    user_options = install.user_options + [
        (
            "debug-gdb",
            None,
            "Add debugging flags for use with GDB (note: reduces performance)",
        ),
        ("num-threads=", "j", "Max number of threads to use for compiling extensions")
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.debug_gdb = 0
        self.num_threads = 0

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        print(self.__dict__)

        import numpy

        extensions = []
        HTSLIB_LIBRARY = ['hts', 'z']
        HTSLIB_LIB_DIRS = [os.environ.get("HTSLIB_LIBRARY_DIR", '/usr/local/lib')]
        HTSLIB_INC_DIRS = [os.environ.get("HTSLIB_INCLUDE_DIR", '/usr/local/include')]

        compile_args = ['-fopenmp']
        if self.debug_gdb:
            compile_args += ["-O0", "-g"]
        else:
            compile_args += ['-O3']

        linker_args = ['-lgomp']

        # if not self.no_lld:
        #     lld_linker = shutil.which('ld.lld')
        #     if lld_linker:
        #         print("Using LLD Linker in place of GNU LD")
        #         os.environ['LDSHARED'] = lld_linker
        #         linker_args += ['-Wl', '--threads' '-Wl', '--thread-count', self.num_threads]

        if sys.platform == 'darwin':
            # os.environ['CLANG_DEFAULT_CXX_STDLIB'] = 'libc++'
            pass
            # compile_args.append('-stdlib=libc++')
        #     linker_args = ['-L/usr/local/opt/llvm/lib'] + linker_args
        else:
            compile_args.append('-std=c++11')

        MAJIQ_INC_DIRS = ['rna_majiq/src/internals']
        MAJIQ_INC_STATS_DIR = ['rna_majiq/src/internals/stats']
        MAJIQ_LIB_DIRS = ['rna_majiq/src/internals']


        NPY_INC_DIRS = [numpy.get_include()]

        # We are currently using a deprecated API (since numpy v1.7)
        # silence deprecation warnings for now
        SILENCE_DEPRECATION = [("NPY_NO_DEPRECATED_API", 0)]

        extensions += [Extension('rna_majiq.src.build', ['rna_majiq/src/build.pyx',
                                                     'rna_majiq/src/internals/io_bam.cpp',
                                                     'rna_majiq/src/internals/grimoire.cpp'],
                                 include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS + HTSLIB_INC_DIRS,
                                 library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                                 libraries=HTSLIB_LIBRARY,
                                 runtime_library_dirs=HTSLIB_LIB_DIRS + MAJIQ_LIB_DIRS,
                                 extra_compile_args=compile_args,
                                 extra_link_args=linker_args,
                                 define_macros=SILENCE_DEPRECATION,
                                 language='c++')]

        extensions += [Extension('rna_majiq.src.calc_psi', ['rna_majiq/src/calc_psi.pyx', 'rna_majiq/src/internals/psi.cpp'],
                                 include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                                 library_dirs=MAJIQ_LIB_DIRS,
                                 runtime_library_dirs=MAJIQ_LIB_DIRS,
                                 extra_compile_args=compile_args,
                                 extra_link_args=linker_args,
                                 define_macros=SILENCE_DEPRECATION,
                                 language='c++')]

        extensions += [Extension('rna_majiq.src.deltapsi', ['rna_majiq/src/deltapsi.pyx', 'rna_majiq/src/internals/psi.cpp'],
                                 include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS,
                                 library_dirs=MAJIQ_LIB_DIRS,
                                 runtime_library_dirs=MAJIQ_LIB_DIRS,
                                 extra_compile_args=compile_args,
                                 extra_link_args=linker_args,
                                 define_macros=SILENCE_DEPRECATION,
                                 language='c++')]

        extensions += [Extension('rna_majiq.src.indpnt', ['rna_majiq/src/indpnt.pyx', 'rna_majiq/src/internals/psi.cpp'],
                                 include_dirs=MAJIQ_INC_DIRS + NPY_INC_DIRS + MAJIQ_INC_STATS_DIR,
                                 library_dirs=MAJIQ_LIB_DIRS,
                                 runtime_library_dirs=MAJIQ_LIB_DIRS,
                                 extra_compile_args=compile_args,
                                 extra_link_args=linker_args,
                                 define_macros=SILENCE_DEPRECATION,
                                 language='c++')]

        extensions += [Extension('rna_majiq.src.io', ['rna_majiq/src/io.pyx'], include_dirs=NPY_INC_DIRS + MAJIQ_INC_DIRS,
                                 define_macros=SILENCE_DEPRECATION,
                                 extra_compile_args=compile_args, extra_link_args=linker_args, language='c++')]

        extensions += [Extension('rna_majiq.c.splice_graph_sql',
                                 ['rna_majiq/c/splice_graph_sql.pyx', 'rna_majiq/c/sqlite3.c'],
                                 language='c++', include_dirs=NPY_INC_DIRS, extra_compile_args=compile_args,
                                 extra_link_args=linker_args)]

        from Cython.Build import cythonize

        self.distribution.ext_modules = cythonize(extensions, language_level=3, nthreads=int(self.num_threads))

        install.run(self)



setup(
    name = 'rna_majiq',
    version = VERSION,
    description="MAJIQ",
    author='BioCiphers Lab',
    author_email='majiq@biociphers.org',
    url='https://biociphers.org',
    keywords=['rna', 'splicing', 'psi', 'splicegraph'],
    license='No License',
    packages=find_packages(),
    setup_requires = setup_requirements(),
    install_requires = requirements(),
    include_package_data = True,
    zip_safe = False,
    entry_points = {'console_scripts': ['majiq = rna_majiq.run_majiq:main']},
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

