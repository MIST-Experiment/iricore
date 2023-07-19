``iricore`` documentation
=========================
The ``iricore`` is the Python wrapper for the
`International Reference Ionosphere <https://irimodel.org/>`_. Currently, the ``iricore`` provides access to
IRI-2016 (last updated: 10/13/2021) and IRI-2020 (last updated: 05/10/2023) versions.
The functionality of the package is limited since only the ``OUTF`` output array from the IRI is processed by
the ``iricore``. The use of the ``OARR`` array is not implemented, which, apart from additional output, disables
the possibility of direct user input for IRI parameters.

.. note::
    The ``iricore`` package is not actively maintained. However, the package can be updated on reasonable
    requests, such as:

    * Updating source IRI version after major releases;
    * Fixing critical bugs that affect all users;
    * Updating the Python version.

    For such inquiries, please contact vadym.bidula@gmail.com

Contents
--------

.. toctree::
    :maxdepth: 2
    :glob:

    installation
    quickstart
    user_guide/index
    reference
