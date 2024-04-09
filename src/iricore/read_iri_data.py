import os
from datetime import datetime

import fortranformat as ff
import numpy as np

_iri_cfd = os.path.dirname(os.path.abspath(__file__))
_datadir = os.path.join(_iri_cfd, f"data/")


def read_apf107():
    with open(os.path.join(_datadir, "index/apf107.dat"), "r") as file:
        lines = file.readlines()
    nlines = len(lines)
    linereader = ff.FortranRecordReader("(13I3, 3F5.1)")
    data = np.array([linereader.read(line_) for line_ in lines], order='F', dtype=np.float32)

    f107d = data[:, 13]
    f107_81 = data[:, 14]
    f107_365 = data[:, 15]
    af107 = np.zeros((27000, 3), order='F', dtype=np.float32)
    af107[:nlines, 0] = f107d
    af107[:nlines, 1] = np.where(f107_81 < -4, f107d, f107_81)
    af107[:nlines, 2] = np.where(f107_365 < -4, f107d, f107_365)
    aap = np.zeros((27000, 9), order='F', dtype=np.int32)
    aap[:nlines, :] = data[:, 3:12].astype(np.int32)
    ly, lm, ld = data[nlines - 1, :3].astype(np.int8)
    last_date = datetime(ly + 2000, lm, ld)
    return (aap, af107, nlines), last_date
