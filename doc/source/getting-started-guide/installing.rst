.. _installing:

Installation
============

Full installation instructions are available after accepting appropriate
license for your use case. For academic use, please see the `majiq download page`_. For commercial use,
you will need to contact us, please refer to `majiq commercial`_.


Required dependencies and tech stack
---------------------

Majiq and Voila are primarily python based tools. In general, the best support for running both the analysis portion
and the visualizer comes from using a Linux-based operating system, however, usage has also been successfully tested
using OSX and windows. In general, before starting, you should have python3.8_ installed, and make sure that both
**python** and
the included **pip** package installation tool are functioning.

In addition, Majiq has one third party library dependency: htslib_ ; Using version **1.10** or higher is highly recommended.
If you are able to administrate the system, it is
simple to install with package managers on Linux-like systems:

- debian-based
    - # apt install libhts-dev
- rpm-based
    - # yum install htslib-devel

Otherwise, the library itself is relatively painless to download and compile manually, either for system-wise or standard
user use. It has no other dependencies at all. Download htslib_ and then follow the standard C compiling steps:

- $ tar xf <downloaded htslib>
- $ cd <untar'd folder>
- $ ./configure --prefix=/htslib/install/location
- $ make
- $ make install

Then, before running any other majiq installation steps, export the following two environment variables (if htslib was
installed to a non-default / non-system-wide location):

- $ export HTSLIB_LIBRARY_DIR=/htslib/install/location/lib
- $ export HTSLIB_INCLUDE_DIR=/htslib/install/location/include
- $ <remaining majiq installation commands>



Usage with Conda
--------

Majiq installation will also work with the 'anaconda' environment manager. (more information about usage or limitations
with conda here)

.. _python3.8: https://www.python.org/downloads/release/python-380/
.. _htslib: http://www.htslib.org/download/
.. _majiq download page: https://majiq.biociphers.org/app_download/
.. _majiq commercial: https://majiq.biociphers.org/commercial.php