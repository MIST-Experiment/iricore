import unittest
from datetime import datetime

import numpy as np

import src.iricore as iricore


# TODO: Fix 2020 tests
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

    def test_iri_ed16(self):
        ed16diff = (self.target_16.edens - self.ref_16.edens) / self.ref_16.edens
        assert np.all(np.where(np.isnan(ed16diff), 0, ed16diff) < 1e-2)

    def test_iri_ed20(self):
        ed20diff = (self.target_20.edens - self.ref_20.edens) / self.ref_20.edens
        assert np.all(np.where(np.isnan(ed20diff), 0, ed20diff) < 1e-2)

    def test_iri_temp16(self):
        nt16diff = (self.target_16.ntemp - self.ref_16.ntemp) / self.ref_16.ntemp
        it16diff = (self.target_16.itemp - self.ref_16.itemp) / self.ref_16.itemp
        et16diff = (self.target_16.etemp - self.ref_16.etemp) / self.ref_16.etemp
        assert np.all(np.where(np.isnan(nt16diff), 0, nt16diff) < 1e-2)
        assert np.all(np.where(np.isnan(it16diff), 0, it16diff) < 1e-2)
        assert np.all(np.where(np.isnan(et16diff), 0, et16diff) < 1e-2)

    def test_iri_temp20(self):
        nt20diff = (self.target_20.ntemp - self.ref_20.ntemp) / self.ref_20.ntemp
        it20diff = (self.target_20.itemp - self.ref_20.itemp) / self.ref_20.itemp
        et20diff = (self.target_20.etemp - self.ref_20.etemp) / self.ref_20.etemp
        assert np.all(np.where(np.isnan(nt20diff), 0, nt20diff) < 1e-2)
        assert np.all(np.where(np.isnan(it20diff), 0, it20diff) < 1e-2)
        assert np.all(np.where(np.isnan(et20diff), 0, et20diff) < 1e-2)

    def test_iri_ions16(self):
        o16diff = (self.target_16.o - self.ref_16.o)
        h16diff = (self.target_16.h - self.ref_16.h)
        he16diff = (self.target_16.he - self.ref_16.he)
        o216diff = (self.target_16.o2 - self.ref_16.o2)
        no16diff = (self.target_16.no - self.ref_16.no)
        cluster16diff = (self.target_16.cluster - self.ref_16.cluster)
        n16diff = (self.target_16.n - self.ref_16.n)
        assert np.all(np.where(np.isnan(o16diff), 0, o16diff) < 1e-1)
        assert np.all(np.where(np.isnan(h16diff), 0, h16diff) < 1e-1)
        assert np.all(np.where(np.isnan(he16diff), 0, he16diff) < 1e-1)
        assert np.all(np.where(np.isnan(o216diff), 0, o216diff) < 1e-1)
        assert np.all(np.where(np.isnan(no16diff), 0, no16diff) < 1e-1)
        assert np.all(np.where(np.isnan(cluster16diff), 0, cluster16diff) < 1e-1)
        assert np.all(np.where(np.isnan(n16diff), 0, n16diff) < 1e-1)

    def test_iri_ions20(self):
        o20diff = (self.target_20.o - self.ref_20.o)
        h20diff = (self.target_20.h - self.ref_20.h)
        he20diff = (self.target_20.he - self.ref_20.he)
        o220diff = (self.target_20.o2 - self.ref_20.o2)
        no20diff = (self.target_20.no - self.ref_20.no)
        n20diff = (self.target_20.n - self.ref_20.n)

        assert np.all(np.where(np.isnan(o20diff), 0, o20diff) < 1e-1)
        assert np.all(np.where(np.isnan(h20diff), 0, h20diff) < 1e-1)
        assert np.all(np.where(np.isnan(he20diff), 0, he20diff) < 1e-1)
        assert np.all(np.where(np.isnan(o220diff), 0, o220diff) < 1e-1)
        assert np.all(np.where(np.isnan(no20diff), 0, no20diff) < 1e-1)
        assert np.all(np.where(np.isnan(n20diff), 0, n20diff) < 1e-1)
