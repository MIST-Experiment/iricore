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
from .iri import iri2016, iri2020, _iri_cfd, iri, indices_uptodate
from .iri_flags import get_jf
from .modules.ion_tools import srange, R_EARTH
from .raytracing import raytrace
from .read_iri_data import read_apf107

_APF107_DATA, _LAST_DATE = read_apf107()

def _clean_ne_for_tec(ne):
    if np.any(np.isnan(ne)) or np.any(np.isinf(ne)):
        warnings.warn("NANs or INFs found in the IRI output; setting them to 0.", stacklevel=2)
    ne = np.where(np.isnan(ne), 0, ne)
    ne = np.where(np.isinf(ne), 0, ne)
    ne_ = ne - np.roll(ne, -1)
    ne = np.where(ne_ > 1e3, np.roll(ne, 1), ne)
    return ne


def _integrate_ne(ne: np.ndarray, hstep: float | np.ndarray) -> np.ndarray:
    # Inefficient for float hstep, but supports arrays
    tec_res = np.sum(ne * hstep, axis=-1) * 1e3 * 1e-16
    if np.any(tec_res < 0) or np.any(tec_res > 200):
        warnings.warn("Unusual TEC value found.", stacklevel=2)
    return tec_res


def _calc_h_steps(height: np.ndarray, el: float) -> np.ndarray:
    step = np.empty(height.shape)
    step[1:-1] = (height[2:] - height[:-2]) / 2
    step[0] = height[1] - height[0]
    step[-1] = height[-1] - height[-2]
    return step / np.cos(np.deg2rad(90 - el))


def _slant2steps(rslant: np.ndarray):
    step = np.empty(rslant.shape)
    step[1:-1] = (rslant[2:] - rslant[:-2]) / 2
    step[0] = rslant[1] - rslant[0]
    step[-1] = rslant[-1] - rslant[-2]
    return step


def vtec(dt: datetime, lat: float | np.ndarray, lon: float | np.ndarray, hbot: float = 90,
         htop: float = 2000, hstep: float = 0.5,
         version: Literal[16, 20] = DEFAULT_IRI_VERSION,
         jf: np.ndarray | str = None, **kwargs):
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

    :return: Vertical TEC.
    """
    indices_uptodate(dt)
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
        iri_res = iri(dt, stage_ranges[i], lat, lon, version, jf, **kwargs)
        tec += _integrate_ne(_clean_ne_for_tec(iri_res.edens), hstep)
    return tec


def stec(el: float, az: float, dt: datetime, lat: float, lon: float, hobs: float = 0, hbot: float = 90,
         htop: float = 2000, npoints: int = 1000, heights: np.ndarray = None,
         version: Literal[16, 20] = DEFAULT_IRI_VERSION,
         jf: np.ndarray | str = None, return_details: bool = False, **kwargs) -> float | Sequence:
    """
    Slant TEC calculated by integrating the electron density on the line of sight. This function
    is not a part of the source IRI code.

    :param el: elevation of observation in [deg].
    :param az: azimuth of observation in [deg].
    :param dt: time of observation.
    :param lat: Geographical latitude.
    :param lon: Geographical longitude.
    :param hobs: Height of the observer above the sea level.
    :param hbot: Bottom height limit for integration in [km].
    :param htop: Upper height limit for integration in [km].
    :param npoints: Number of points to integrate.
    :param heights: Overrides the height grid with a custom array of height in [km]
    :param version: IRI version number.
    :param jf: Array of JF parameters or string for predefined JF arrays to be used in the IRI_SUB function. See
               :func:`iricore.get_jf` for details. If not specified otherwise, the default IRI JF array will be used.
    :param return_details: If True - also returns history with ray position, corresponding electron density and the oarr
                           array. Dict keys: ['lat', 'lon', 'h', 'edens', 'oarr']
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

    :return: Slant TEC.
    """
    indices_uptodate(dt)
    if hbot < 60 or htop > 2000:
        raise ValueError("The limits of integration cannot exceed (60, 2000) km.")
    el = np.array(el)
    az = np.array(az)

    if heights is None:
        # hstep = (htop - hbot) / npoints
        heights = np.linspace(hbot, htop, npoints)
    # Calculating input parameters (assuming Earth=sphere)
    rslant = srange(np.deg2rad(90 - el), heights * 1e3)
    ell = pm.Ellipsoid(R_EARTH, R_EARTH)
    slat, slon, _ = pm.aer2geodetic(az, el, rslant, lat, lon, hobs, ell=ell)
    h_hist = np.empty((*el.shape, heights.size))
    h_hist[..., :] = heights
    history = {'lat': slat, 'lon': slon, 'h': h_hist, 'ds': rslant * 1e-3}
    history['ds'][1:] = history['ds'][1:] - history['ds'][:-1]

    if not isinstance(jf, np.ndarray):
        jf = jf or "default_edens"
        if isinstance(jf, str):
            jf = get_jf(jf)
    else:
        if isinstance(jf, np.ndarray) and jf.size != 50:
            raise ValueError("Length of jf array must be 50")

    iri_res, oarr = _call_stec(dt, heights, slat, slon, jf, version, **kwargs)
    ne = iri_res[0].transpose()
    ne = ne.reshape((len(heights), -1))[:, 0]
    history['edens'] = ne
    history['oarr'] = oarr

    # Data postprocessing (fixing non-physical IRI output)
    ne = _clean_ne_for_tec(ne)
    tec_res = _integrate_ne(ne, _slant2steps(rslant * 1e-3))
    if return_details:
        return tec_res, history
    return tec_res


def _call_stec(dt: datetime, heights: Sequence[float], lat: Sequence[float], lon: Sequence[float],
               jf: np.ndarray, version: Literal[16, 20] = DEFAULT_IRI_VERSION, **kwargs):
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

    # Handling manual user input
    if len(kwargs) > 0:
        if "oarr0" in kwargs.keys():
            jf[7] = 0
            oarr[0] = kwargs["oarr0"]
        if "oarr1" in kwargs.keys():
            jf[8] = 0
            oarr[1] = kwargs["oarr1"]
        if "oarr14" in kwargs.keys():
            jf[9] = 0
            oarr[14] = kwargs["oarr14"][0]
            oarr[15] = kwargs["oarr14"][1]
        if "oarr2" in kwargs.keys():
            jf[12] = 0
            oarr[2] = kwargs["oarr2"]
        if "oarr3" in kwargs.keys():
            jf[13] = 0
            oarr[3] = kwargs["oarr3"]
        if "oarr4" in kwargs.keys():
            jf[14] = 0
            oarr[4] = kwargs["oarr4"]
        if "oarr5" in kwargs.keys():
            jf[15] = 0
            oarr[5] = kwargs["oarr5"]
        if "oarr32" in kwargs.keys():
            jf[16] = 0
            oarr[32] = kwargs["oarr32"]
        if "oarr40" in kwargs.keys():
            jf[24] = 0
            oarr[40] = kwargs["oarr40"]
        if "oarr38" in kwargs.keys():
            jf[26] = 0
            oarr[38] = kwargs["oarr38"]
        if "oarr45" in kwargs.keys():
            jf[31] = 0
            oarr[45] = kwargs["oarr45"]
        if "oarr9" in kwargs.keys():
            jf[42] = 0
            oarr[9] = kwargs["oarr9"]
        if "oarr34" in kwargs.keys():
            jf[43] = 0
            oarr[34] = kwargs["oarr34"]

    datadir = jpath(_iri_cfd, 'data')
    datadir_bytes = bytes(datadir, 'utf-8')

    aap, af107, nlines = _APF107_DATA

    iricore.stec_(as_ctypes(jf), byref(c_bool(jmag)), f_lat, f_lon, f_heights, byref(hsize), byref(iyyyy), byref(mmdd),
                  byref(dhour), as_ctypes(oarr),
                  iri_res.ctypes.data_as(POINTER(c_float)), datadir_bytes, byref(c_int(len(datadir))),
                  aap.ctypes.data_as(POINTER(c_float)), af107.ctypes.data_as(POINTER(c_float)), byref(c_int(nlines)))

    iri_res = np.ascontiguousarray(iri_res)
    return iri_res, oarr


def refstec(el: float, az: float, dt: datetime, lat: float, lon: float, freq: float,
            hobs: float = 0, hbot: float = 90,
            htop: float = 2000, npoints: int = 1000, heights: np.ndarray = None,
            jf: np.ndarray | str = None, return_hist: bool = False) -> float | Sequence:
    """
    Slant TEC calculated by integrating the electron density on the line of sight. This function
    is not a part of the source IRI code.

    :param el: elevation of observation in [deg].
    :param az: azimuth of observation in [deg].
    :param dt: time of observation.
    :param lat: Geographical latitude.
    :param lon: Geographical longitude.
    :param freq: Frequency of observation.
    :param hobs: Height of the observer above the sea level.
    :param hbot: Bottom height limit for integration in [km].
    :param htop: Upper height limit for integration in [km].
    :param npoints: Number of points to integrate.
    :param heights: Overrides the height grid with a custom array of height in [km]
    :param jf: Array of JF parameters or string for predefined JF arrays to be used in the IRI_SUB function. See
               :func:`iricore.get_jf` for details. If not specified otherwise, the default IRI JF array will be used.
    :param return_hist: If True - also returns history with ray position and corresponding electron density.
                        Dict keys: ['lat', 'lon', 'h', 'edens']
    :return: Slant TEC.
    """
    indices_uptodate(dt)
    if hbot < 60 or htop > 2000:
        raise ValueError("The limits of integration cannot exceed (60, 2000) km.")

    if heights is None:
        # hstep = (htop - hbot) / npoints
        heights = np.linspace(hbot, htop, npoints)
    # Calculating input parameters (assuming Earth=sphere)
    rslant = srange(np.deg2rad(90 - el), heights * 1e3)
    ell = pm.Ellipsoid(R_EARTH, R_EARTH)
    slat, slon, _ = pm.aer2geodetic(az, el, rslant, lat, lon, hobs, ell=ell)

    if not isinstance(jf, np.ndarray):
        jf = jf or "default_edens"
        if isinstance(jf, str):
            jf = get_jf(jf)
    else:
        if isinstance(jf, np.ndarray) and jf.size != 50:
            raise ValueError("Length of jf array must be 50")
    # TODO: pass jf to raytrace
    pos = (lat, lon, hobs)
    history = raytrace(el, az, freq, pos, heights, dt)
    ne = history["edens"]

    # Data postprocessing (fixing non-physical IRI output)
    ne = _clean_ne_for_tec(ne)
    steps = np.empty(history['ds'].shape)
    steps[1:-1] = (history['ds'][1:-1] + history['ds'][2:]) / 2
    steps[0] = history['ds'][1]
    steps[-1] = history['ds'][-1]

    tec_res = _integrate_ne(ne, steps)
    if return_hist:
        return tec_res, history
    return tec_res
