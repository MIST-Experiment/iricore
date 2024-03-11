# A simple 2D raytracing. Code copied from dionpy
import datetime
from typing import Tuple

import numpy as np
import pymap3d as pm
from pymap3d import Ellipsoid

from .iri import iri
from .modules.ion_tools import srange, refr_index, refr_angle, plasfreq, R_EARTH

_ROUND_ELL = Ellipsoid(R_EARTH, R_EARTH)
_LIGHT_SPEED = 2.99792458e8  # in [m/s]


def _raytrace_sublayer(lat_ray, lon_ray, h_ray, h_next, alt_cur, az, freq, d_theta, ref_ind, n_sublayer, nlayers, dt,
                       theta_ref=None):
    # Distance from current position to next layer
    r_slant = srange(np.deg2rad(90 - alt_cur), h_next - h_ray, re=R_EARTH + h_ray)
    lat_next, lon_next, _ = pm.aer2geodetic(az, alt_cur, r_slant, lat_ray, lon_ray, h_ray, ell=_ROUND_ELL)
    # The sides of the 1st triangle
    d_cur = R_EARTH + h_ray  # Distance from Earth center to current point
    d_next = R_EARTH + h_next  # Distance from Earth center to layer

    if theta_ref is None:
        # The inclination angle at the interface using law of cosines [rad]
        costheta_inc = (r_slant ** 2 + d_next ** 2 - d_cur ** 2) / (2 * r_slant * d_next)
        assert (costheta_inc <= 1).all(), (f"Cosine of inclination angle cannot be >= 1. Something is wrong with "
                                           f"coordinates at heights {h_ray * 1e-3:.1f}-{h_next * 1e-3:.1f} [km].")
        theta_inc = np.arccos(costheta_inc)
    else:
        # Angle between d_cur and r_slant
        int_angle_rad = np.pi - theta_ref
        # The inclination angle at the i-th interface using law of sines [rad]
        theta_inc = np.arcsin(np.sin(int_angle_rad) * d_cur / d_next)

    # Get IRI info of point
    ed = iri(dt, [h_next / 1e3, h_next / 1e3, 1], lat_next, lon_next).edens
    ed = np.where(ed < 0, 0, ed)

    # Refraction index of the surface
    if n_sublayer == nlayers - 1:
        ref_ind_next = np.ones(alt_cur.shape)
        nan_theta_mask = np.zeros(alt_cur.shape)
    else:
        ref_ind_next = refr_index(ed, freq)
        nan_theta_mask = plasfreq(ed, angular=False) > freq

    # The outgoing angle at the 1st interface using Snell's law
    theta_ref = refr_angle(ref_ind, ref_ind_next, theta_inc)
    inf_theta_mask = np.abs((ref_ind / ref_ind_next * np.sin(theta_inc))) > 1
    d_theta += theta_ref - theta_inc
    alt_next = np.rad2deg(np.pi / 2 - theta_ref)
    return lat_next, lon_next, h_next, d_theta, alt_next, ref_ind_next, theta_ref, ed, nan_theta_mask, inf_theta_mask, r_slant


def raytrace(
        alt: float | np.ndarray,
        az: float | np.ndarray,
        freq: float | np.ndarray,
        pos: Tuple[float, float, float],
        heights: np.ndarray,
        dt: datetime.datetime
) -> dict:
    # Initialization of variables
    freq *= 1e6
    alt_cur = np.array(alt)
    az = np.array(az)
    heights = heights * 1e3  # in [m]
    nlayers = heights.size

    delta_theta = 0 * alt_cur
    delta_theta_hist = np.empty((*alt_cur.shape, nlayers))
    inf_theta_mask = 0 * alt_cur
    nan_theta_mask = 0 * alt_cur

    history = {
        'edens': np.zeros((*alt_cur.shape, nlayers)),
        'lat': np.zeros((*alt_cur.shape, nlayers)),
        'lon': np.zeros((*alt_cur.shape, nlayers)),
        'h': np.zeros((*alt_cur.shape, nlayers)),
        'ds': np.zeros((*alt_cur.shape, nlayers))
    }

    # Init values for the first sub-layer
    ref_ind_cur = np.ones(alt_cur.shape)
    lat_ray, lon_ray, h_ray = pos
    theta_ref = None
    for i in range(nlayers):
        # Tracing change in position due to refraction
        (lat_ray, lon_ray, h_ray, delta_theta, alt_cur, ref_ind_cur, theta_ref, ed,
         nt_mask, it_mask, ds) = _raytrace_sublayer(
            lat_ray, lon_ray, h_ray, heights[i], alt_cur, az, freq, delta_theta, ref_ind_cur, i, nlayers, dt, theta_ref)
        delta_theta_hist[..., i] = delta_theta
        nan_theta_mask += nt_mask
        inf_theta_mask += it_mask
        history['edens'][..., i] = ed
        history['lat'][..., i] = lat_ray
        history['lon'][..., i] = lon_ray
        history['h'][..., i] = h_ray / 1e3
        history['ds'][..., i] = ds / 1e3
        history['lat'][..., i] = np.where(inf_theta_mask == 0, history['lat'][..., i], np.inf)
        history['lon'][..., i] = np.where(inf_theta_mask == 0, history['lon'][..., i], np.inf)
        history['lat'][..., i] = np.where(nan_theta_mask == 0, history['lat'][..., i], np.nan)
        history['lon'][..., i] = np.where(nan_theta_mask == 0, history['lon'][..., i], np.nan)

        if inf_theta_mask:
            history['lat'][..., i:] = np.inf
            history['lon'][..., i:] = np.inf
            history['edens'][..., i] = np.nan
            history['h'][..., i] = np.nan
            history['ds'][..., i] = np.nan
            return history

    return history


def raytrace_star(args):
    """
    For parallel calculations
    """
    return raytrace(*args)
