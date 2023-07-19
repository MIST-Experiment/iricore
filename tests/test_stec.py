import unittest
import numpy as np
from datetime import datetime
import src.iricore as iricore

class TestSTEC(unittest.TestCase):
    def test_stec():
        dt = datetime(2022, 10, 5, 12)
        LAT, LON = 79.3798, -90.99885
        el = np.linspace(0, 90, 10)
        az = np.full(10, 45)
        tecs = np.empty(10)
        test_tecs = np.load("test_tecs.npy")

        for i in range(el.size):
            tecs[i] = iricore.stec(
                el.ravel()[i],
                az.ravel()[i],
                dt,
                (LAT, LON, 0),
            )
        assert np.all(np.isclose(tecs, test_tecs))
