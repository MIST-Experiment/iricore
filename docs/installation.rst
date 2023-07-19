Installation
============


During the installation, the source IRI code will be compiled automatically in the background, which requires
CMAKE and any Fortran compiler pre-installed.

Linux
-----

* Install a Fortran compiler, e.g. `GFortran <https://gcc.gnu.org/wiki/GFortran>`_ and CMAKE:
    .. code-block::

        sudo apt install gfortran cmake

* Use ``pip`` or any other Python package manager to install the ``iricore``:

    .. code-block::

        python3 -m pip install iricore

Windows
-------
Although the Windows platform is not supported by ``iricore``, it is not impossible to use it on Windows.
If you are on Windows, consider installing `WSL <https://docs.microsoft.com/en-us/windows/wsl/install>`_ and then follow
the steps for Linux.
If that is not an option, then:

* Gather all of your luck;
* Install a Fortran compiler. For example, you can follow `this guide <https://masuday.github.io/fortran_tutorial/install_gfortran_windows.html>`_;
* Install `CMAKE <https://cmake.org/>`_;

* Use ``pip`` or any other Python package manager to (attempt to) install the ``iricore``:

    .. code-block::

        python3 -m pip install iricore

* Fix the errors that appear during installation if there are any.
