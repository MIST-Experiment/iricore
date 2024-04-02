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
    :param oarr: Additional output data. First index axis corresponds to different coordinates, for example,
                 oarr[0] will correspond to lat[0] and lon[0] and so on. The second index axis has size 100 and maps
                 the extra output information. If only one coordinate was calculated, the first axis is discarded.
                        * oarr[0] = NMF2/M-3
                        * oarr[1] = HMF2/KM
                        * oarr[2] = NMF1/M-3
                        * oarr[3] = HMF1/KM
                        * oarr[4] = NME/M-3
                        * oarr[5] = HME/KM
                        * oarr[6] = NMD/M-3
                        * oarr[7] = HMD/KM
                        * oarr[8] = HHALF/KM
                        * oarr[9] = B0/KM
                        * oarr[10] = VALLEY-BASE/M-3
                        * oarr[11] = VALLEY-TOP/KM
                        * oarr[12] = TE-PEAK/K
                        * oarr[13] = TE-PEAK HEIGHT/KM
                        * oarr[14] = TE-MOD(300KM)
                        * oarr[15] = TE-MOD(400KM)/K
                        * oarr[16] = TE-MOD(600KM)
                        * oarr[17] = TE-MOD(1400KM)/K
                        * oarr[18] = TE-MOD(3000KM)
                        * oarr[19] = TE(120KM)=TN=TI/K
                        * oarr[20] = TI-MOD(430KM)
                        * oarr[21] = X/KM, WHERE TE=TI
                        * oarr[22] = SOL ZENITH ANG/DEG
                        * oarr[23] = SUN DECLINATION/DEG
                        * oarr[24] = DIP/deg
                        * oarr[25] = DIP LATITUDE/deg
                        * oarr[26] = MODIFIED DIP LAT.
                        * oarr[27] = Geographic latitude
                        * oarr[28] = sunrise/dec. hours
                        * oarr[29] = sunset/dec. hours
                        * oarr[30] = ISEASON (1=spring)
                        * oarr[31] = Geographic longitude
                        * oarr[32] = Rz12
                        * oarr[33] = Covington Index
                        * oarr[34] = B1
                        * oarr[35] = M(3000)F2
                        * oarr[36] = TEC/m-2
                        * oarr[37] = TEC_top/TEC*100.
                        * oarr[38] = gind (IG12)
                        * oarr[39] = F1 probability
                        * oarr[40] = F10.7 daily
                        * oarr[41] = c1 (F1 shape)
                        * oarr[42] = daynr
                        * oarr[43] = equatorial vertical ion drift in m/s
                        * oarr[44] = foF2_storm/foF2_quiet
                        * oarr[45] = F10.7_81
                        * oarr[46] = foE_storm/foE_quiet
                        * oarr[47] = spread-F probability
                        * oarr[48] = Geomag. latitude
                        * oarr[49] = Geomag. longitude
                        * oarr[50] = ap at current time
                        * oarr[51] = daily ap
                        * oarr[52] = invdip/degree
                        * oarr[53] = MLT-Te
                        * oarr[54] = CGM-latitude
                        * oarr[55] = CGM-longitude
                        * oarr[56] = CGM-MLT
                        * oarr[57] = CGM lat eq. aurl bodry
                        * oarr[58] = CGM-lati(MLT=0)
                        * oarr[59] = CGM-lati for MLT=1
                        * oarr[60] = CGM-lati(MLT=2)
                        * oarr[61] = CGM-lati for MLT=3
                        * oarr[62] = CGM-lati(MLT=4)
                        * oarr[63] = CGM-lati for MLT=5
                        * oarr[64] = CGM-lati(MLT=6)
                        * oarr[65] = CGM-lati for MLT=7
                        * oarr[66] = CGM-lati(MLT=8)
                        * oarr[67] = CGM-lati for MLT=9
                        * oarr[68] = CGM-lati(MLT=10)
                        * oarr[69] = CGM-lati for MLT=11
                        * oarr[70] = CGM-lati(MLT=12)
                        * oarr[71] = CGM-lati for MLT=13
                        * oarr[72] = CGM-lati(MLT=14)
                        * oarr[73] = CGM-lati for MLT=15
                        * oarr[74] = CGM-lati(MLT=16)
                        * oarr[75] = CGM-lati for MLT=17
                        * oarr[76] = CGM-lati(MLT=18)
                        * oarr[77] = CGM-lati for MLT=19
                        * oarr[78] = CGM-lati(MLT=20)
                        * oarr[79] = CGM-lati for MLT=21
                        * oarr[80] = CGM-lati(MLT=22)
                        * oarr[81] = CGM-lati for MLT=23
                        * oarr[82] = Kp at current time
                        * oarr[83] = magnetic declination
                        * oarr[84] = L-value
                        * oarr[85] = dipole moment
                        * oarr[86] = SAX300
                        * oarr[87] = SUX300
                        * oarr[88] = HNEA
                        * oarr[89] = HNEE
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
    oarr: np.ndarray

    @staticmethod
    def _extract_iri_row(row, arr_shape):
        ncoord, nalt = arr_shape
        row = row.transpose().reshape((ncoord, -1))[:, :nalt]
        return np.squeeze(np.where(row < 0, np.nan, row))

    @classmethod
    def from_raw(cls, outf_out: np.ndarray, oarr_out: np.ndarray, lat: np.ndarray, lon: np.ndarray,
                 height_range: [float, float, float]):
        height = np.arange(height_range[0], height_range[1] + height_range[2] / 2, height_range[2])
        arr_shape = (lat.size, height.size)
        return cls(
            lat=lat,
            lon=lon,
            height=height,
            edens=cls._extract_iri_row(outf_out[0], arr_shape),
            ntemp=cls._extract_iri_row(outf_out[1], arr_shape),
            itemp=cls._extract_iri_row(outf_out[2], arr_shape),
            etemp=cls._extract_iri_row(outf_out[3], arr_shape),
            o=cls._extract_iri_row(outf_out[4], arr_shape),
            h=cls._extract_iri_row(outf_out[5], arr_shape),
            he=cls._extract_iri_row(outf_out[6], arr_shape),
            o2=cls._extract_iri_row(outf_out[7], arr_shape),
            no=cls._extract_iri_row(outf_out[8], arr_shape),
            cluster=cls._extract_iri_row(outf_out[9], arr_shape),
            n=cls._extract_iri_row(outf_out[10], arr_shape),
            firi_output=cls._extract_iri_row(outf_out[13], arr_shape),
            oarr=np.squeeze(oarr_out.swapaxes(0, 1))
        )
