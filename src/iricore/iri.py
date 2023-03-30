from __future__ import annotations

import glob
import os
from ctypes import *
from datetime import datetime
from os.path import join as jpath
from typing import Sequence

import numpy as np
from numpy.ctypeslib import as_ctypes

from .read_iri_data import readapf107

_iri_cfd = os.path.dirname(os.path.abspath(__file__))


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

IRI_VERSIONS = ['16', '20']


def _call_iri_sub(dt: datetime, alt_range: [float, float, float], lats: Sequence[float], lons: Sequence[float],
                  jf: np.ndarray, version: int = 20, aap: np.ndarray | None = None,
                  af107: np.ndarray | None = None, nlines: int | None = None):
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

    if aap is None and af107 is None and nlines is None:
        aap, af107, nlines = readapf107(version)

    # ==================================================================================================================
    iricore.iricore_(as_ctypes(jf), byref(c_bool(jmag)), glat, glon, byref(gsize), byref(iyyyy), byref(mmdd),
                     byref(dhour), byref(heibeg), byref(heiend), byref(heistp), as_ctypes(oarr),
                     iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))),
                     aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines)))
    # ==================================================================================================================

    iri_res = np.ascontiguousarray(iri_res)
    return iri_res


def _call_stec(dt: datetime, heights: Sequence[float], lat: Sequence[float], lon: Sequence[float],
               jf: np.ndarray, version: int = 16, aap: np.ndarray | None = None,
               af107: np.ndarray | None = None, nlines: int | None = None):
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

    if aap is None and af107 is None and nlines is None:
        aap, af107, nlines = readapf107(version)

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
        replace_missing: float = np.nan, version: int = 20, aap: np.ndarray | None = None,
        af107: np.ndarray | None = None, nlines: int | None = None) -> dict:
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[2, 3, 4, 5, 11, 20, 21, 22, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_iri_sub(dt, alt_range, lats, lons, jf, version, aap, af107, nlines)
    ne = _extract_data(iri_res, 0, len(lats), alt_range, replace_missing)
    te = _extract_data(iri_res, 3, len(lats), alt_range, replace_missing)
    return {'ne': ne, 'te': te}


def IRI_etemp_only(dt: datetime, alt_range: [float, float, float], lats: Sequence[float], lons: Sequence[float],
                   replace_missing: float = np.nan, version: int = 20, aap: np.ndarray | None = None,
                   af107: np.ndarray | None = None, nlines: int | None = None) -> dict:
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[0, 2, 3, 4, 5, 11, 20, 21, 22, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_iri_sub(dt, alt_range, lats, lons, jf, version, aap, af107, nlines)
    te = _extract_data(iri_res, 3, len(lats), alt_range, replace_missing)
    return {'te': te}


def STEC(dt: datetime, heights: Sequence[float], lats: Sequence[float], lons: Sequence[float],
         version: int = 20, aap: np.ndarray | None = None,
         af107: np.ndarray | None = None, nlines: int | None = None) -> np.ndarray:
    jf = np.ones(50, dtype=np.int32, order="F")
    jf[[1, 2, 3, 4, 5, 11, 20, 21, 22, 25, 27, 28, 29, 33, 34, 35, 36, 46]] = 0
    iri_res = _call_stec(dt, heights, lats, lons, jf, version, aap, af107, nlines)
    ne = iri_res[0].transpose()
    ne = ne.reshape((len(heights), -1))[:, 0]
    return ne
