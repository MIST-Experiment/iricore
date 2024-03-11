import unittest
from datetime import datetime

import numpy as np

import src.iricore as iricore


class TestIRI(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestIRI, self).__init__(*args, **kwargs)
        dt = datetime(2020, 4, 11, 10)
        lat = 10
        lon = 110
        altrange = [0, 1000, 10]
        alt = np.linspace(altrange[0], altrange[1] + 1, altrange[2])
        self.target_16 = iricore.iri(dt, altrange, lat, lon, version=16)
        self.target_20 = iricore.iri(dt, altrange, lat, lon, version=20)

        ref_data_16 = np.loadtxt("data/test_iri2016_data.txt", unpack=True)
        ref_data_16 = np.where(ref_data_16 == -1, np.nan, ref_data_16)
        self.ref_16 = iricore.IRIOutput(
            lat, lon, alt,
            ref_data_16[1] * 1e6,
            ref_data_16[3],
            ref_data_16[4],
            ref_data_16[5],
            ref_data_16[6] / 10,
            ref_data_16[8] / 10,
            ref_data_16[9] / 10,
            ref_data_16[10] / 10,
            ref_data_16[11] / 10,
            ref_data_16[12] / 10,
            ref_data_16[7] / 10,
            None,
        )

        ref_data_20 = np.load("data/iri20_data.npy", )
        ref_data_20 = np.where(ref_data_20 == -1, np.nan, ref_data_20)
        self.ref_20 = iricore.IRIOutput(
            lat, lon, alt,
            ref_data_20[0],
            ref_data_20[1],
            ref_data_20[2],
            ref_data_20[3],
            ref_data_20[4],
            ref_data_20[5],
            ref_data_20[6],
            ref_data_20[7],
            ref_data_20[8],
            None,
            ref_data_20[9],
            None,
        )

    # TODO: Fix tests

    def test_iri_ed16(self):
        np.testing.assert_allclose(self.target_16.edens, self.ref_16.edens, rtol=1e-2)

    def test_iri_ed20(self):
        np.testing.assert_allclose(self.target_20.edens, self.ref_20.edens, rtol=1e-2)

    def test_iri_temp16(self):
        np.testing.assert_allclose(self.target_16.ntemp, self.ref_16.ntemp, rtol=1e-2)
        np.testing.assert_allclose(self.target_16.itemp, self.ref_16.itemp, rtol=1e-2)
        np.testing.assert_allclose(self.target_16.etemp, self.ref_16.etemp, rtol=1e-2)

    def test_iri_temp20(self):
        np.testing.assert_allclose(self.target_20.ntemp, self.ref_20.ntemp, rtol=1e-2)
        np.testing.assert_allclose(self.target_20.itemp, self.ref_20.itemp, rtol=1e-2)
        np.testing.assert_allclose(self.target_20.etemp, self.ref_20.etemp, rtol=1e-2)
    #
    # def test_iri_ions16(self):
    #     # print(self.target_16.o)
    #     # print(self.ref_16.o)
    #     np.testing.assert_allclose(self.target_16.o, self.ref_16.o, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.h, self.ref_16.h, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.he, self.ref_16.he, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.o2, self.ref_16.o2, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.no, self.ref_16.no, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.cluster, self.ref_16.cluster, atol=1e-1)
    #     np.testing.assert_allclose(self.target_16.n, self.ref_16.n, atol=1e-1)
    #
    #
    # def test_iri_ions20(self):
    #     print(self.target_20.o)
    #     print(self.ref_20.o)
    #     np.testing.assert_allclose(self.target_20.o, self.ref_20.o, atol=1e-1)
    #     np.testing.assert_allclose(self.target_20.h, self.ref_20.h, atol=1e-1)
    #     np.testing.assert_allclose(self.target_20.he, self.ref_20.he, atol=1e-1)
    #     np.testing.assert_allclose(self.target_20.o2, self.ref_20.o2, atol=1e-1)
    #     np.testing.assert_allclose(self.target_20.no, self.ref_20.no, atol=1e-1)
    #     np.testing.assert_allclose(self.target_20.n, self.ref_20.n, atol=1e-1)
