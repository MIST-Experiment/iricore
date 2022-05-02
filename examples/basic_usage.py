import iricore
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

# iricore.IRI() is the main function in the package, which provides access to all the necessary values. First of all,
# it depends on the single datetime object, which defines date and time of the observation
dt = datetime(year=2018, month=2, day=13, hour=6, minute=20)

# The second argument is a range of heights for calculation: [start point in km, end point in km, step in km]. If you
# want to calculate single height, just type the same values for start and end points, e.g. [500, 500, 1]
heights = [100, 500, 10]

# The last two arguments are arrays of latitudes and longitudes. Important! Values from those arrays are matched index
# by index, which means if you specify array of, lets say, 10 latitudes and 10 longitudes the result will be calculated
# for 10 coordinate points - (lat[0], lon[0]), (lat[1], lon[1]), (lat[2], lon[2]) and so on.
lats = np.linspace(0, 90, 100)
lons = np.repeat(0, 100)

result = iricore.IRI(dt, heights, lats, lons)

# Output of iricore.IRI is a dictionary with two keys - 'ne' (electron density) and 'te' (electron temperature). Both
# these keys provide access to 2D numpy array. First index indicates coordinate point, e.g. result['ne'][5] gives
# access to electron density at all specified heights for the coordinate point (lat[5], lon[5]). Second index gives
# access to the heights. Examples:
ind = 5
# I increase height[1] by 1 to include endpoint
height_array = np.arange(heights[0], heights[1]+1, heights[2])
plt.plot(height_array, result['ne'][ind])
plt.xlabel("Height, km")
plt.ylabel(r"$n_e$, $m^{-3}$")
plt.title(r"$n_e$" + f" profile for lat = {lats[ind]:.2f} deg, lon = {lons[ind]:.2f} deg")
plt.show()


ind = 10
plt.plot(lats, result['ne'][:, ind])
plt.xlabel("Latitude, deg")
plt.ylabel(r"$n_e$, $m^{-3}$")
plt.title(r"$n_e$" + f" latitude profile for h = {height_array[ind]:.2f} km")
plt.show()

