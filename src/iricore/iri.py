from __future__ import annotations

import glob
import os
import warnings
from ctypes import *
from datetime import datetime
from os.path import join as jpath
from typing import Sequence, Literal
from .config import IRI_VERSIONS, DEFAULT_VERSION

import numpy as np
import pymap3d as pm
from numpy.ctypeslib import as_ctypes

from .modules import srange, R_EARTH
from .read_iri_data import readapf107

_iri_cfd = os.path.dirname(os.path.abspath(__file__))

IRI_DATA = readapf107(20)

def _import_libs():
    iri2016 = np.ctypeslib.load_library("libiri2016", _iri_cfd)
    iri2020 = np.ctypeslib.load_library("libiri2020", _iri_cfd)
    return iri2016, iri2020


def _move_libs():
    parent_dir = os.path.abspath(jpath(_iri_cfd, '..'))
    libs = glob.glob(jpath(parent_dir, 'libiri*'))
    for lib in libs:
        parent, file = os.path.split(lib)
        os.rename(lib, jpath(parent, 'iricore', file))


try:
    iri2016, iri2020 = _import_libs()
except OSError:
    _move_libs()
    try:
        iri2016, iri2020 = _import_libs()
    except OSError:
        raise ImportError("Could not import IRI libraries. Please make sure you have installed the package correctly.")


def _call_iri_sub(dt: datetime, alt_range: [float, float, float], lats: Sequence[float], lons: Sequence[float],
                  jf: np.ndarray, version: Literal[16, 20] = DEFAULT_VERSION):
    if version == 16:
        iricore = iri2016
    elif version == 20:
        iricore = iri2020
    else:
        raise ValueError(f"Available IRI versions: " + ", ".join(IRI_VERSIONS))

    lats = np.asarray(lats)
    lons = np.asarray(lons)
    if not len(lats) == len(lons):
        raise ValueError("Lengths of latitude and longitude arrays must be equal.")

    jmag = False
    iyyyy = c_int(dt.year)
    mmdd = c_int(100 * dt.month + dt.day)
    dhour = c_float(dt.hour + dt.minute / 60 + dt.second / 3600 + 25.)
    glat = as_ctypes(np.array(lats, dtype=np.float32, order="F"))
    glon = as_ctypes(np.array(lons, dtype=np.float32, order="F"))
    gsize = c_int(len(lats))
    heibeg = c_float(alt_range[0])
    heiend = c_float(alt_range[1])
    heistp = c_float(alt_range[2])
    oarr = np.zeros(100, dtype=np.float32, order="F")
    iri_res = np.zeros((20, 1000, len(lats)), dtype=np.float32, order='F')

    datadir = jpath(_iri_cfd, 'data/data' + str(version))
    datadir_bytes = bytes(datadir, 'utf-8')
    aap, af107, nlines = IRI_DATA

    # ==================================================================================================================
    iricore.iricore_(as_ctypes(jf), byref(c_bool(jmag)), glat, glon, byref(gsize), byref(iyyyy), byref(mmdd),
                     byref(dhour), byref(heibeg), byref(heiend), byref(heistp), as_ctypes(oarr),
                     iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))),
                     aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines)))
    # ==================================================================================================================

    iri_res = np.ascontiguousarray(iri_res)
    return iri_res


def _call_stec(dt: datetime, heights: Sequence[float], lat: Sequence[float], lon: Sequence[float],
               jf: np.ndarray, version: Literal[16, 20] = DEFAULT_VERSION):
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

    datadir = jpath(_iri_cfd, 'data/data' + str(version))
    datadir_bytes = bytes(datadir, 'utf-8')

    aap, af107, nlines = IRI_DATA

    # ==================================================================================================================
    iricore.stec_(as_ctypes(jf), byref(c_bool(jmag)), f_lat, f_lon, f_heights, byref(hsize), byref(iyyyy), byref(mmdd),
                  byref(dhour), as_ctypes(oarr),
                  iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))),
                  aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines)))
    # ==================================================================================================================

    iri_res = np.ascontiguousarray(iri_res)
    return iri_res


def _extract_data(iri_res: np.ndarray, index: int, ncoord: int, alt_range: [float, float, float],
                  replace_missing: float = np.nan):
    nalts = int((alt_range[1] - alt_range[0]) / alt_range[2]) + 1
    res = iri_res[index].transpose()
    res = res.reshape((ncoord, -1))[:, :nalts]
    return np.where(res < 0, replace_missing, res)


def IRI(dt: datetime, alt_range: [float, float, float], lats: Sequence[float], lons: Sequence[float],
        replace_missing: float = np.nan, version: Literal[16, 20] = DEFAULT_VERSION) -> dict:
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[2, 3, 4, 5, 11, 20, 21, 22, 23, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_iri_sub(dt, alt_range, lats, lons, jf, version)
    ne = _extract_data(iri_res, 0, len(lats), alt_range, replace_missing)
    te = _extract_data(iri_res, 3, len(lats), alt_range, replace_missing)
    return {'ne': ne, 'te': te}


def IRI_etemp_only(dt: datetime, alt_range: [float, float, float], lats: Sequence[float], lons: Sequence[float],
                   replace_missing: float = np.nan, version: Literal[16, 20] = DEFAULT_VERSION) -> dict:
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[0, 2, 3, 4, 5, 11, 20, 21, 22, 23, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_iri_sub(dt, alt_range, lats, lons, jf, version)
    te = _extract_data(iri_res, 3, len(lats), alt_range, replace_missing)
    return {'te': te}


def stec(alt: float, az: float, dt: datetime, position: Sequence[float, float, float], hbot: float = 90,
         htop: float = 2000, npoints: int = 500,
         version: Literal[16, 20] = DEFAULT_VERSION, debug: bool = False) -> float | Sequence:
    """
    :param alt: altitude (elevation) of observation in [deg].
    :param az: azimuth of observation in [deg].
    :param dt: time of observation.
    :param position: sequence containing geographical latilude, longitude and altitude
                     above sea level in [deg, deg, m].
    :param hbot: Bottom height limit for integration in [km].
    :param htop: Upper height limit for integration in [km].
    :param npoints: Number of points to integrate.
    :param version: IRI version number.
    :param debug: If True - also returns IRI Ne output and Ne after post-processing.
    :return: Slant TEC.
    """
    # Calculating input parameters (assuming Earth=sphere)
    hstep = (htop - hbot) / npoints
    heights = np.linspace(hbot, htop, npoints)
    rslant = srange(np.deg2rad(90 - alt), heights*1e3)
    ell = pm.Ellipsoid(R_EARTH, R_EARTH)
    slat, slon, _ = pm.aer2geodetic(az, alt, rslant, *position, ell=ell)

    # IRI call and data extraction
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[1, 2, 3, 4, 5, 11, 20, 21, 22, 23, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_stec(dt, heights, slat, slon, jf, version)
    ne = iri_res[0].transpose()
    ne = ne.reshape((len(heights), -1))[:, 0]

    # Data postprocessing (fixing non-physical IRI output)
    ne_debug = ne.copy()
    if np.any(np.isnan(ne)) or np.any(np.isinf(ne)):
        warnings.warn("NANs or INFs found in the IRI output; setting them to 0.", stacklevel=2)
    ne = np.where(np.isnan(ne), 0, ne)
    ne = np.where(np.isinf(ne), 0, ne)
    ne_ = ne - np.roll(ne, -1)
    ne = np.where(ne_ > 1e3, np.roll(ne, 1), ne)
    stec_res = np.sum(ne) * hstep * 1e3 * 1e-16
    if stec_res < 0 or stec_res > 100:
        warnings.warn("Unusual sTEC value found.", stacklevel=2)
    if debug:
        return [stec_res, ne_debug, ne]
    return stec_res
