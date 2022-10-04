# This script compares time of execution and results from the iricore and original iri2016_old package. To run this script
# you will need two additional packages - iri2016_old and tqdm. You can install them with
# >>> python -m pip install tqdm iri2016_old
# During the first run the iri2016_old will compile itself, which can take some time.

import iricore
import iri2016
from datetime import datetime
import numpy as np
from time import time
from tqdm import tqdm

dt = datetime(year=2018, month=1, day=1, hour=1)    # Picking random date

N = 1000    # Number of iterations

t_start = time()
iri16_res = []
# Iterating over array of latitudes, keeping longitude the same
for lat in tqdm(np.linspace(0, 50, N)):
    iri16_res.append(iri2016.IRI(dt, [500, 510, 1], lat, 0))

t_iri16 = time() - t_start
print(f"Time spent with iri2016: {t_iri16:.3f} s")

t_start = time()
# In case of iricore you don't need to set up an explicit loop - passing array of latitudes and longitudes is enough
irico_res = iricore.IRI(dt, [500, 510, 1], np.linspace(0, 50, N), np.repeat(0, N))
t_irico = time() - t_start
print(f"Time spent with iricore: {t_irico:.3f} s")

# Making sure that both packages provide the same result
iri_2016_ne = np.array([res.ne.data for res in iri16_res])
iri_2016_te = np.array([res.Te.data for res in iri16_res])
print("Same n_e?\t", np.isclose(iri_2016_ne, irico_res['ne']).all())
print("Same T_e?\t", np.isclose(iri_2016_te, irico_res['te']).all())
print(f"Performance boost ~ {int(t_iri16/t_irico)}x")
