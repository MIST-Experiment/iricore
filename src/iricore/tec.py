from __future__ import annotations

import warnings
from ctypes import *
from datetime import datetime
from os.path import join as jpath
from typing import Sequence, Literal

import numpy as np
import pymap3d as pm
from numpy.ctypeslib import as_ctypes

from .config import IRI_VERSIONS, DEFAULT_IRI_VERSION
from .iri import iri2016, iri2020, _iri_cfd, iri
from .iri_flags import get_jf
from .modules import srange, R_EARTH


def _clean_ne_for_tec(ne):
    if np.any(np.isnan(ne)) or np.any(np.isinf(ne)):
        warnings.warn("NANs or INFs found in the IRI output; setting them to 0.", stacklevel=2)
    ne = np.where(np.isnan(ne), 0, ne)
    ne = np.where(np.isinf(ne), 0, ne)
    ne_ = ne - np.roll(ne, -1)
    ne = np.where(ne_ > 1e3, np.roll(ne, 1), ne)
    return ne


def _integrate_ne(ne, hstep):
    tec_res = np.sum(ne, axis=-1) * hstep * 1e3 * 1e-16
    if np.any(tec_res < 0) or np.any(tec_res > 100):
        warnings.warn("Unusual TEC value found.", stacklevel=2)
    return tec_res


def vtec(dt: datetime, lat: float | np.ndarray, lon: float | np.ndarray, hbot: float = 90,
         htop: float = 2000, hstep: float = 0.5,
         version: Literal[16, 20] = DEFAULT_IRI_VERSION,
         jf: np.ndarray | str = None):
    """
    Vertical TEC calculated by integrating the electron density on the line of sight. This function
    is not a part of the source IRI code.

    :param dt: time of observation.
    :param lat: Geographical latitude.
    :param lon: Geographical longitude.
    :param hbot: Bottom height limit for integration in [km].
    :param htop: Upper height limit for integration in [km].
    :param hstep: Height step for integration in [km].
    :param version: IRI version number.
    :param jf: Array of JF parameters or string for predefined JF arrays to be used in the IRI_SUB function. See
               :func:`iricore.get_jf` for details. If not specified otherwise, the default IRI JF array will be used.
    :return: Slant TEC.
    """

    if hbot < 60 or htop > 2000:
        raise ValueError("The limits of integration cannot exceed (60, 2000) km.")

    if not isinstance(jf, np.ndarray):
        jf = jf or "default_edens"
        if isinstance(jf, str):
            jf = get_jf(jf)
    else:
        if isinstance(jf, np.ndarray) and jf.size != 50:
            raise ValueError("Length of jf array must be 50")

    npoints = int(np.ceil((htop - hbot) / hstep))
    nstages = int(np.ceil(npoints / 1000))
    stage_ranges = [[
        hbot + hstep * (1000 * i),
        hbot + hstep * (1000 * (i + 1) - 1),
        hstep
    ] for i in range(nstages)]
    stage_ranges[-1][1] = htop

    tec = 0.
    for i in range(nstages):
        iri_res = iri(dt, stage_ranges[i], lat, lon, version, jf)
        tec += _integrate_ne(_clean_ne_for_tec(iri_res.edens), hstep)
    return tec


def stec(el: float, az: float, dt: datetime, lat: float, lon: float, height: float = 0, hbot: float = 90,
         htop: float = 2000, npoints: int = 1000,
         version: Literal[16, 20] = DEFAULT_IRI_VERSION,
         jf: np.ndarray | str = None, _return_ne: bool = False) -> float | Sequence:
    """
    Slant TEC calculated by integrating the electron density on the line of sight. This function
    is not a part of the source IRI code.

    :param el: elevation of observation in [deg].
    :param az: azimuth of observation in [deg].
    :param dt: time of observation.
    :param lat: Geographical latitude.
    :param lon: Geographical longitude.
    :param height: Height of the observer above the sea level.
    :param hbot: Bottom height limit for integration in [km].
    :param htop: Upper height limit for integration in [km].
    :param npoints: Number of points to integrate.
    :param version: IRI version number.
    :param jf: Array of JF parameters or string for predefined JF arrays to be used in the IRI_SUB function. See
               :func:`iricore.get_jf` for details. If not specified otherwise, the default IRI JF array will be used.
    :param _return_ne: If True - also returns IRI Ne output and Ne after and before
                       post-processing (used for debugging).
    :return: Slant TEC.
    """
    if hbot < 60 or htop > 2000:
        raise ValueError("The limits of integration cannot exceed (60, 2000) km.")

    # Calculating input parameters (assuming Earth=sphere)
    hstep = (htop - hbot) / npoints
    heights = np.linspace(hbot, htop, npoints)
    rslant = srange(np.deg2rad(90 - el), heights * 1e3)
    ell = pm.Ellipsoid(R_EARTH, R_EARTH)
    slat, slon, _ = pm.aer2geodetic(az, el, rslant, lat, lon, height, ell=ell)

    if not isinstance(jf, np.ndarray):
        jf = jf or "default_edens"
        if isinstance(jf, str):
            jf = get_jf(jf)
    else:
        if isinstance(jf, np.ndarray) and jf.size != 50:
            raise ValueError("Length of jf array must be 50")

    iri_res = _call_stec(dt, heights, slat, slon, jf, version)
    ne = iri_res[0].transpose()
    ne = ne.reshape((len(heights), -1))[:, 0]

    # Data postprocessing (fixing non-physical IRI output)
    ne_debug = ne.copy()
    ne = _clean_ne_for_tec(ne)
    tec_res = _integrate_ne(ne, hstep)
    if _return_ne:
        return tec_res, ne, ne_debug
    return tec_res


def _call_stec(dt: datetime, heights: Sequence[float], lat: Sequence[float], lon: Sequence[float],
               jf: np.ndarray, version: Literal[16, 20] = DEFAULT_IRI_VERSION):
    if version == 16:
        iricore = iri2016
    elif version == 20:
        iricore = iri2020
    else:
        raise ValueError(f"Available IRI versions: " + ", ".join(IRI_VERSIONS))

    heights = np.asarray(heights)
    jmag = False
    iyyyy = c_int(dt.year)
    mmdd = c_int(100 * dt.month + dt.day)
    dhour = c_float(dt.hour + dt.minute / 60 + dt.second / 3600 + 25.)
    f_heights = as_ctypes(np.array(heights, dtype=np.float32, order="F"))
    f_lat = as_ctypes(np.array(lat, dtype=np.float32, order="F"))
    f_lon = as_ctypes(np.array(lon, dtype=np.float32, order="F"))
    hsize = c_int(len(heights))
    oarr = np.zeros(100, dtype=np.float32, order="F")
    iri_res = np.zeros((20, 1000, len(heights)), dtype=np.float32, order='F')

    datadir = jpath(_iri_cfd, 'data')
    datadir_bytes = bytes(datadir, 'utf-8')

    # aap, af107, nlines = IRI_DATA

    iricore.stec_(as_ctypes(jf), byref(c_bool(jmag)), f_lat, f_lon, f_heights, byref(hsize), byref(iyyyy), byref(mmdd),
                  byref(dhour), as_ctypes(oarr),
                  iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))))

    iri_res = np.ascontiguousarray(iri_res)
    return iri_res
