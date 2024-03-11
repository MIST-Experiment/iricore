from __future__ import annotations

import numpy as np

np.seterr(invalid="ignore")

R_EARTH = 6378100.0  # in [m]


def srange(
        theta: float | np.ndarray, alt: float | np.ndarray, re: float = 6378100
) -> float | np.ndarray:
    """
    :param theta: Zenith angle in [rad].
    :param alt: Altitude in [m].
    :param re: Radius of the Earth in [m].
    :return: Distance in meters from the telescope to the point (theta, alt) in [m]
    """
    r = -re * np.cos(theta) + np.sqrt(
        (re * np.cos(theta)) ** 2 + alt ** 2 + 2 * alt * re
    )
    return r


def plasfreq(n_e: float | np.ndarray, angular: bool = True) -> float | np.ndarray:
    """
    Angular (omega) plasma frequency.

    :param n_e: Electron density in [m^-3].
    :param angular: If True - return angular frequency, otherwise frequency.
    :return: Plasma frequency of cold electrons in Hz.
    """
    e = 1.60217662e-19
    m_e = 9.10938356e-31
    epsilon0 = 8.85418782e-12
    if np.min(n_e) < 0:
        raise ValueError(
            "Number density cannot be < 0. Most probably iricore does not include data for the specified date. Please "
            "update the library by calling iricore.update()."
        )
    result = np.sqrt((n_e * e ** 2) / (m_e * epsilon0))
    return result if angular else 0.5 * result / np.pi


def refr_index(n_e: float | np.ndarray, freq: float):
    """

    :param n_e: Electron density in [m^-3].
    :param freq: Observational frequency in [Hz].
    :return: Refractive index of the ionosphere from electron density.
    """
    nu_p = plasfreq(n_e, angular=False)
    return np.sqrt(1 - (nu_p / freq) ** 2)


def refr_angle(
        n1: float | np.ndarray,
        n2: float | np.ndarray,
        phi: float | np.ndarray,
) -> float | np.ndarray:
    """
    Snell's law.

    :param n1: Refractive index in previous medium.
    :param n2: Refractive index in current medium.
    :param phi: Angle of incident ray in [rad].
    :return: Outcoming angle in [rad].
    """
    return np.arcsin(n1 / n2 * np.sin(phi))


def trop_refr(el: float | np.ndarray, h: float) -> float | np.ndarray:
    """
    Approximation of the refraction in the troposphere recommended by the ITU-R:
    https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.834-9-201712-I!!PDF-E.pdf

    :param el: Elevation angle in [deg].
    :param h: Height of instrument in [km].
    :return: Change of the angle theta due to tropospheric refraction (in radians).
    """
    t1 = 1.314 + 0.6437 * el + 0.02869 * el ** 2
    t2 = 0.2305 + 0.09428 * el + 0.01096 * el ** 2
    return 1 / (t1 + h * t2 + h ** 2 * 0.08583)
