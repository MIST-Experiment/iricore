import os

import fortranformat as ff
import numpy as np

_iri_cfd = os.path.dirname(os.path.abspath(__file__))


def readapf107():
    datadir = os.path.join(_iri_cfd, f"data/")
    with open(os.path.join(datadir, "index/apf107.dat"), "r") as file:
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
    af107[:nlines, 1] = np.where(f107_365 < -4, f107d, f107_365)

    aap = np.zeros((27000, 9), order='F', dtype=np.float32)
    aap[:nlines, :] = data[:, 3:12].astype(np.int32)
    return aap, af107, nlines
