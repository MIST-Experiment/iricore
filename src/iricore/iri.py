from __future__ import annotations

import os
from ctypes import *
from datetime import datetime
from os.path import join as jpath
from typing import Sequence, Literal, Annotated

import numpy as np
from numpy.ctypeslib import as_ctypes

from .config import IRI_VERSIONS, DEFAULT_IRI_VERSION
from .iri_flags import get_jf
from .irioutput import IRIOutput

_iri_cfd = os.path.dirname(os.path.abspath(__file__))

# TODO: Fix data reading from Python
# TODO: Check if data works for all files
# IRI_DATA = readapf107()

try:
    iri2016 = np.ctypeslib.load_library("libiri2016", _iri_cfd)
    iri2020 = np.ctypeslib.load_library("libiri2020", _iri_cfd)
except OSError:
    raise ImportError(
        """
        Could not import IRI libraries, most probably they were not compiled during installation. You can try 
        compiling them manually:
        1) Locate the package installation directory and open it in terminal
        2) Run: 
            >>> cmake .
            >>> make
        """
    )


def iri(dt: datetime, altrange: Annotated[Sequence[float], 3], lat: float | Sequence[float],
        lon: float | Sequence[float], version: Literal[16, 20] = DEFAULT_IRI_VERSION,
        jf: np.ndarray | str = None) -> IRIOutput:
    """
    The main function of the package representing the IRI_SUB routine from the IRI source code. Provides access
    to all IRI calculations through the wide range of parameters.

    :param dt: Date and time - must be a single value.
    :param altrange: Range of altitudes in the form (alt_start, alt_stop, alt_step), all in [km].
    :param lat: Geographical latitude.
    :param lon: Geographical longitude.
    :param version: Version of the IRI to use:

                    * ``version=16`` - IRI-2016
                    * ``version=20`` - IRI-2020

    :param jf: Array of JF parameters or string for predefined JF arrays to be used in the IRI_SUB function. See
               :func:`iricore.get_jf` for details. If not specified otherwise, the default IRI JF array will be used.
    :return: :class:`IRIOutput` dataclass.
    """
    lat = np.atleast_1d(np.asarray(lat))
    lon = np.atleast_1d(np.asarray(lon))
    if not len(lat) == len(lon):
        raise ValueError("Lengths of latitude and longitude arrays must be equal.")

    if not isinstance(jf, np.ndarray):
        jf = jf or "default"
        if isinstance(jf, str):
            jf = get_jf(jf)
    else:
        if isinstance(jf, np.ndarray) and jf.size != 50:
            raise ValueError("Length of jf array must be 50")

    if version == 16:
        iricore = iri2016
    elif version == 20:
        iricore = iri2020
    else:
        raise ValueError(f"Available IRI versions: " + ", ".join(IRI_VERSIONS))

    if (altrange[1] - altrange[0]) / altrange[2] + 1 > 1000:
        raise ValueError("The specified altitude range and step require more that 1000 points for calculation, "
                         "which exceeds IRI limitations. Please consider breaking calculation into chunks.")

    # Converting parameters to Fortran types
    jmag = False
    iyyyy = c_int(dt.year)
    mmdd = c_int(100 * dt.month + dt.day)
    dhour = c_float(dt.hour + dt.minute / 60 + dt.second / 3600 + 25.)
    glat = as_ctypes(np.array(lat, dtype=np.float32, order="F"))
    glon = as_ctypes(np.array(lon, dtype=np.float32, order="F"))
    gsize = c_int(len(lat))
    heibeg = c_float(altrange[0])
    heiend = c_float(altrange[1])
    heistp = c_float(altrange[2])
    oarr = np.zeros(100, dtype=np.float32, order="F")
    iri_res = np.zeros((20, 1000, len(lat)), dtype=np.float32, order='F')

    datadir = jpath(_iri_cfd, 'data')
    datadir_bytes = bytes(datadir, 'utf-8')
    # aap, af107, nlines = IRI_DATA

    # Calling IRI_SUB
    iricore.iricore_(as_ctypes(jf), byref(c_bool(jmag)), glat, glon, byref(gsize), byref(iyyyy), byref(mmdd),
                     byref(dhour), byref(heibeg), byref(heiend), byref(heistp), as_ctypes(oarr),
                     iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))),
                     # aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines))
                     )

    iri_res = np.ascontiguousarray(iri_res)
    return IRIOutput.from_raw(iri_res, lat, lon, altrange)
