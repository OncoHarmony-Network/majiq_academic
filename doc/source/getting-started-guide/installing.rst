.. _installing:

Installation
============

Full installation instructions are available after accepting appropriate
license for your use case.
For academic use, please see the `majiq download page`_.
For commercial use, you will need to contact us, please refer to
`majiq commercial`_.


Required dependencies
---------------------

MAJIQ and VOILA are best supported with a Linux-based operating system.
However, they have also successfully been tested/used these tools on MacOS and
Windows.

MAJIQ requires for installation:

- Python (3.8 or later)
- setuptools (45 or later)
- HTSlib_ (1.10 or later)
- C++ compiler with C++11 support (e.g. gcc >= 4.8.1)

In general, before starting, you should have python3.8_ installed, with both
**python** and the included **pip** package installation tool.


HTSlib installation
~~~~~~~~~~~~~~~~~~~

MAJIQ's depends on HTSlib_ (1.10 or later) in order to efficiently parse BAM
files from RNA-seq experiments.

If you are have admin rights, you can install system-wide with package managers
on Linux:

.. role:: bash(code)
   :language: bash

- debian-based: :bash:`apt install libhts-dev`,
- rpm-based: :bash:`yum install htslib-devel`.


Otherwise, you can install HTSlib_ from source.
For example, to install htslib-1.13 to ``~/install/htslib-1.13``, you could run:

.. code-block:: bash

   # download htslib 1.13 archive to current directory
   curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
   tar -xf htslib-1.13.tar.bz2  # extract archive
   cd htslib-1.13  # change into htslib source directory
   # configure, make, and install htslib to ~/install/htslib-1.13
   ./configure --prefix=$HOME/install/htslib-1.13
   make
   make install


See ``INSTALL`` from the HTSlib_ sources for detailed instructions.


MAJIQ installation
~~~~~~~~~~~~~~~~~~

If HTSlib was not installed to the default system-wide location, MAJIQ
installers need to know where it was installed.
This is done by setting the environment variables ``HTSLIB_LIBRARY_DIR`` and
``HTSLIB_INCLUDE_DIR``.

For example, with previous instructions, from the same directory as the README:

.. code-block:: bash

   # change to where library/include directories are installed
   export HTSLIB_LIBRARY_DIR=$HOME/install/htslib-1.13/lib
   export HTSLIB_INCLUDE_DIR=$HOME/install/htslib-1.13/include

Afterwards, install MAJIQ following the full instructions available after
accepting the license agreement
(e.g., :bash:`pip install <url to MAJIQ repository>`).


.. _python3.8: https://www.python.org/downloads/release/python-380/
.. _HTSlib: http://www.htslib.org/download/
.. _majiq download page: https://majiq.biociphers.org/app_download/
.. _majiq commercial: https://majiq.biociphers.org/commercial.php
