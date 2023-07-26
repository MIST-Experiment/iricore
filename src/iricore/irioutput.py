from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class IRIOutput:
    """
    Dataclass containing :func:`iricore.iri` output. If lat/lon are 1D arrays, then the
    IRI returned values will a 2D array with shape (lat.size, alt.size). If lat/lon
    are a single number, the returned values will be an 1D array with shape (alt.size).

    :param lat: Array of latitudes corresponding to the IRI output (0 axis).
    :param lon: Array of longitudes corresponding to the IRI output (0 axis).
    :param height: Array of heights corresponding to the IRI output (1 axis).
    :param edens: Electron density in [m-3].
    :param ntemp: Neutral temperature in [K].
    :param itemp: Ion temperature in [K].
    :param etemp: Electron temperature in [K].
    :param o: O+ ion density in [%](default) or [m-3].
    :param h: H+ ion density in [%](default) or [m-3].
    :param he: He+ ion density in [%](default) or [m-3].
    :param o2: O2+ ion density in [%](default) or [m-3].
    :param no: NO+ ion density in [%](default) or [m-3].
    :param cluster: Cluster ion density in [%](default) or [m-3].
    :param n: N+ ion density in [%](default) or [m-3].
    :param firi_output: If FIRI was manually selected for the D-region, the ``firi_output`` will be an 1D array
                        containing:

                        * [0:10] standard IRI-Ne for 60,65,..,110km;
                        * [11:21] Friedrich (FIRI) model at these heights;
                        * [22:32] standard Danilov (SW=0, WA=0);
                        * [33:43] for minor Stratospheric Warming (SW=0.5);
                        * [55:65] weak Winter Anomaly (WA=0.5) conditions;
                        * [44:54] for major Stratospheric Warming (SW=1);
                        * [66:76] strong Winter Anomaly (WA=1) conditions.
    """
    lat: np.ndarray
    lon: np.ndarray
    height: np.ndarray
    edens: np.ndarray
    ntemp: np.ndarray
    itemp: np.ndarray
    etemp: np.ndarray
    o: np.ndarray
    h: np.ndarray
    he: np.ndarray
    o2: np.ndarray
    no: np.ndarray
    cluster: np.ndarray
    n: np.ndarray
    firi_output: np.ndarray

    @staticmethod
    def _extract_iri_row(row, arr_shape):
        ncoord, nalt = arr_shape
        row = row.transpose().reshape((ncoord, -1))[:, :nalt]
        return np.squeeze(np.where(row < 0, np.nan, row))

    @classmethod
    def from_raw(cls, output: np.ndarray, lat: np.ndarray, lon: np.ndarray, height_range: [float, float, float]):
        height = np.arange(height_range[0], height_range[1] + height_range[2] / 2, height_range[2])
        arr_shape = (lat.size, height.size)
        return cls(
            lat=lat,
            lon=lon,
            height=height,
            edens=cls._extract_iri_row(output[0], arr_shape),
            ntemp=cls._extract_iri_row(output[1], arr_shape),
            itemp=cls._extract_iri_row(output[2], arr_shape),
            etemp=cls._extract_iri_row(output[3], arr_shape),
            o=cls._extract_iri_row(output[4], arr_shape),
            h=cls._extract_iri_row(output[5], arr_shape),
            he=cls._extract_iri_row(output[6], arr_shape),
            o2=cls._extract_iri_row(output[7], arr_shape),
            no=cls._extract_iri_row(output[8], arr_shape),
            cluster=cls._extract_iri_row(output[9], arr_shape),
            n=cls._extract_iri_row(output[10], arr_shape),
            firi_output=cls._extract_iri_row(output[13], arr_shape),
        )
