# Function to identify grid cells with multiple points from a list
# here, stations with same cells.

import numpy as np
import pandas as pd
import h3
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
import cartopy.crs as ccrs

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'

stations = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1)
atlantic_lon_index = np.where(np.logical_and(stations['Longitude'] >= -100, stations['Longitude'] <= 20))

atlantic_lon = np.take(stations['Longitude'], atlantic_lon_index[0])
atlantic_lat = np.take(stations['Latitude'], atlantic_lon_index[0])

hex_list = [h3.geo_to_h3(y, x, 3) for x, y in zip(atlantic_lon.ravel(), atlantic_lat.ravel())]

print(len(hex_list))

print(len(np.unique(hex_list)))
unique_hex_list, count = np.unique(hex_list, return_counts=True)

dups_index = count > 1
dups_hex = unique_hex_list[dups_index]

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='50m')


def plot_hex(hash, c):
    polygons = h3.h3_set_to_multi_polygon([hash], geo_json=False)
    p = Polygon(polygons[0][0])
    y, x = p.exterior.xy
    ax.plot(x, y, color=c, linewidth=0.5)


for h in unique_hex_list:
    plot_hex(str(h), 'royalblue')

for h in dups_hex:
    plot_hex(str(h), 'red')

ax.scatter(atlantic_lon, atlantic_lat, s=1, c='g')

# for xy in zip(atlantic_lon, atlantic_lat):
#     ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data')
#     print('(%s, %s)' % xy)

plt.show()
