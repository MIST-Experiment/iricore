from __future__ import annotations

import numpy as np

R_EARTH = 6378100.


def srange(theta: float | np.ndarray, alt: float | np.ndarray) -> float | np.ndarray:
    """
    :param theta: Zenith angle in [rad].
    :param alt: Altitude in [m].
    :param re: Radius of the Earth in [m].
    :return: Distance in meters from the telescope to the point (theta, alt)
    """
    if isinstance(theta, np.ndarray) and isinstance(alt, np.ndarray):
        raise ValueError("Only one input parameter can be a numpy array.")
    r = -R_EARTH * np.cos(theta) + np.sqrt(
        (R_EARTH * np.cos(theta)) ** 2 + alt ** 2 + 2 * alt * R_EARTH
    )
    return r

