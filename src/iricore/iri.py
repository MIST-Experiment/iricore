from __future__ import annotations

import os
from ctypes import *
from datetime import datetime
from os.path import join as jpath
from typing import Sequence, Literal, Annotated

import numpy as np
from numpy.ctypeslib import as_ctypes

from .config import IRI_VERSIONS, DEFAULT_IRI_VERSION
from .data_update import update
from .iri_flags import get_jf
from .irioutput import IRIOutput
from .read_iri_data import read_apf107

_iri_cfd = os.path.dirname(os.path.abspath(__file__))

_APF107_DATA, _LAST_DATE = read_apf107()


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


def indices_uptodate(dt: datetime):
    global _APF107_DATA, _LAST_DATE
    if dt > _LAST_DATE:
        print("Requested date is not covered in database. Updating indices.")
        if update():
            _APF107_DATA, _LAST_DATE = read_apf107()
            if dt > _LAST_DATE:
                raise RuntimeError("Requested indices are not available yet. The update latency is usually 2-3 days.")
        else:
            raise RuntimeError("Cannot update indices. Check internet connection.")


def iri(dt: datetime, altrange: Annotated[Sequence[float], 3], lat: float | Sequence[float],
        lon: float | Sequence[float], version: Literal[16, 20] = DEFAULT_IRI_VERSION,
        jf: np.ndarray | str = None, **kwargs) -> IRIOutput:
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
    :param kwargs: Use it to enter user input values for the IRI calculation. Names of parameters are
                   'oarr'+'[oarr index you want to modify]'. Available parameters are:

                   * **oarr0** - float - user input for foF2/MHz or NmF2/m-3
                   * **oarr1** - float - user input for hmF2/km or M(3000)F2
                   * **oarr2** - float - user input for foF1/MHz or NmF1/m-3
                   * **oarr3** - float - user input for hmF1/km
                   * **oarr4** - float - user input for foE/MHz or NmE/m-3
                   * **oarr5** - float - user input for hmE/km
                   * **oarr9** - float - user input for B0
                   * **oarr14** - (float, float) - user input for Ne(300km), Ne(400km)/m-3. Use oarr14[...]=-1 if one of
                     these values is not available. If jf(22)==False then Ne(300km), Ne(550km)/m-3.
                   * **oarr32** - float - user input for Rz12
                   * **oarr34** - float - user input for B1
                   * **oarr38** - float - user input for IG12
                   * **oarr40** - float - user input for daily F10.7 index (make sure to also specify oarr45, otherwise it
                     will be copied from oarr40)
                   * **oarr45** - float - user input for 81-day avg F10.7

    :return: :class:`iricore.IRIOutput` dataclass.
    """
    indices_uptodate(dt)
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
        raise ValueError("The specified altitude range and step require more than 1000 points for calculation, "
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
    oarr_out = np.zeros((100, len(lat)), dtype=np.float32, order="F")
    outf_out = np.zeros((20, 1000, len(lat)), dtype=np.float32, order='F')

    # Handling manual user input
    if len(kwargs) > 0:
        if "oarr0" in kwargs.keys():
            jf[7] = 0
            oarr_out[0, :] = kwargs["oarr0"]
        if "oarr1" in kwargs.keys():
            jf[8] = 0
            oarr_out[1, :] = kwargs["oarr1"]
        if "oarr14" in kwargs.keys():
            jf[9] = 0
            oarr_out[14, :] = kwargs["oarr14"][0]
            oarr_out[15, :] = kwargs["oarr14"][1]
        if "oarr2" in kwargs.keys():
            jf[12] = 0
            oarr_out[2, :] = kwargs["oarr2"]
        if "oarr3" in kwargs.keys():
            jf[13] = 0
            oarr_out[3, :] = kwargs["oarr3"]
        if "oarr4" in kwargs.keys():
            jf[14] = 0
            oarr_out[4, :] = kwargs["oarr4"]
        if "oarr5" in kwargs.keys():
            jf[15] = 0
            oarr_out[5, :] = kwargs["oarr5"]
        if "oarr32" in kwargs.keys():
            jf[16] = 0
            oarr_out[32, :] = kwargs["oarr32"]
        if "oarr40" in kwargs.keys():
            jf[24] = 0
            oarr_out[40, :] = kwargs["oarr40"]
        if "oarr38" in kwargs.keys():
            jf[26] = 0
            oarr_out[38, :] = kwargs["oarr38"]
        if "oarr45" in kwargs.keys():
            jf[31] = 0
            oarr_out[45, :] = kwargs["oarr45"]
        if "oarr9" in kwargs.keys():
            jf[42] = 0
            oarr_out[9, :] = kwargs["oarr9"]
        if "oarr34" in kwargs.keys():
            jf[43] = 0
            oarr_out[34, :] = kwargs["oarr34"]

    datadir = jpath(_iri_cfd, 'data')
    datadir_bytes = bytes(datadir, 'utf-8')
    aap, af107, nlines = _APF107_DATA

    # Calling IRI_SUB
    iricore.iricore_(as_ctypes(jf), byref(c_bool(jmag)), glat, glon, byref(gsize),
                     byref(iyyyy), byref(mmdd), byref(dhour),
                     byref(heibeg), byref(heiend), byref(heistp),
                     # oarr_out.ctypes.data_as(POINTER(c_float)),
                     oarr_out.ctypes.data_as(POINTER(c_float)),
                     outf_out.ctypes.data_as(POINTER(c_float)),
                     datadir_bytes, byref(c_int(len(datadir))),
                     aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines))
                     )

    outf_out = np.ascontiguousarray(outf_out)
    oarr_out = np.ascontiguousarray(oarr_out)
    return IRIOutput.from_raw(outf_out, oarr_out, lat, lon, altrange)
