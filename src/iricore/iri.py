import numpy as np
from typing import Iterable
from datetime import datetime
import os

try:
    from . import iri_fcore as core
except ImportError:
    raise ImportError("Cannot import compiled IRI library.")


class MissingDataError(Exception):
    pass


def IRI(dt: datetime, alt_range: [float, float, float], lats: Iterable[float], lons: Iterable[float],
        replace_missing: float = np.nan) -> dict:
    try:
        from . import iri_fcore as core
    except ImportError:
        raise ImportError("Cannot import compiled IRI library.")

    lats = np.asarray(lats)
    lons = np.asarray(lons)

    if not len(lats) == len(lons):
        raise ValueError("Lengths of latitude and longitude arrays must be equal.")

    jf = np.ones(50, dtype=bool)
    jf[[3, 4, 5, 11, 21, 22, 25, 27, 28, 29, 32, 33, 34, 35, 36]] = 0
    jmag = 0
    mmdd = 100 * dt.month + dt.day
    dhour = dt.hour + dt.minute / 60 + dt.second / 3600
    oarr = np.zeros(100, dtype=np.float32)
    core.mod.iri_res = np.zeros((20, 1000, len(lats)), dtype=np.float32)
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/')
    core.mod.iri_core(jf, jmag, lats, lons, dt.year, mmdd, dhour + 25., alt_range[0], alt_range[1], alt_range[2], oarr,
                      datadir)
    ne = core.mod.iri_res[0].transpose()
    te = core.mod.iri_res[3].transpose()
    nalts = int((alt_range[1] - alt_range[0]) / alt_range[2]) + 1
    ne = ne.reshape((len(lats), -1))[:, :nalts]
    te = te.reshape((len(lats), -1))[:, :nalts]
    res = {
        'ne': np.where(ne < 0, replace_missing, ne),
        'te': np.where(te < 0, replace_missing, te),
    }
    return res
