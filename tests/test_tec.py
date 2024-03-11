import unittest
from datetime import datetime

import numpy as np

import src.iricore as iricore


class TestTEC(unittest.TestCase):
    def test_stec(self):
        dt = datetime(2022, 10, 5, 12)
        LAT, LON = 79.3798, -90.99885
        el = np.linspace(0, 90, 10)
        az = np.full(10, 45)
        tecs = np.empty(10)
        test_tecs = np.load("data/stec_data.npy")

        for i in range(el.size):
            tecs[i] = iricore.stec(
                el.ravel()[i],
                az.ravel()[i],
                dt,
                LAT,
                LON,
                0,
                version=20,
            )
        np.testing.assert_allclose(tecs, test_tecs)

    def test_vtec(self):
        ref_tec = np.load('data/vtec_data.npy')
        dt = datetime(2021, 4, 11, 10)
        lat = np.linspace(0, 90, 10)
        lon = np.linspace(0, 180, 10)
        test_tec = iricore.vtec(dt, lat, lon, version=20)
        np.testing.assert_allclose(ref_tec, test_tec)
