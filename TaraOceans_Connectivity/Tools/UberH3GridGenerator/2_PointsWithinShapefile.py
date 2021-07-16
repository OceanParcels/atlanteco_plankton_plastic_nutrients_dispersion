# Shapefile source: https://www.marineregions.org/gazetteer.php?p=details&id=1902

import shapefile as shp
from shapely.ops import unary_union
from shapely.geometry import shape as shape_geom
from shapely.geometry import Point
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.colors as clr

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
file1 = home_folder + 'iho_Atlantic/iho.shp'
file2 = home_folder + 'iho_Arctic/iho.shp'
# (Doing SO separately as I need to remove particles to the west of 75 W
# file1 = home_folder + 'iho_SO/iho.shp'

# with shp.Reader(file) as sf:
sf1 = shp.Reader(file1)
print(sf1)
sf2 = shp.Reader(file2)
print(sf2)

poly_list = sf1.shapes() + sf2.shapes()

t1 = time.time()
# dissolve multiple polygons to get one outer bound
union = unary_union([shape_geom(s.__geo_interface__) for s in poly_list])  # returns shapely type polygon
print("Time to get the outer Polygon: ", time.time() - t1)

# lat lon from Uber H3
coords = pd.read_csv(home_folder + 'ROMEO_Res6_Centroids.csv')
y = coords['Latitudes']
x = coords['Longitudes']

# exclude points outside the final polygon
# numpy ravel(view) is faster than flatten(copy)
t1 = time.time()
Points = [Point(x, y) for x, y in zip(x.ravel(), y.ravel())]
print("Total number of points: ", len(Points))
Points = unary_union(Points)

# intersection over within or contains- this one return points(MultiPoint)  that lie within the polygon
result = Points.intersection(union)
print("Time to get the points within the outer polygon: ", time.time() - t1)

print("Points within: ", len(result))

final_lats = np.array([a.y for a in result])
final_lons = np.array([a.x for a in result])

# for SO - remove points on the west of 80 deg W
# lons_index = np.where(final_lons >= -75)
# final_lats = final_lats[lons_index]
# final_lons = final_lons[lons_index]

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

fig = plt.figure()
ax = plt.axes()
colormap = clr.ListedColormap(['gainsboro', 'whitesmoke'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

ax.scatter(final_lons, final_lats, s=1, c='g')
plt.show()

coords = {"Latitudes": final_lats, "Longitudes": final_lons}
df = pd.DataFrame(coords, columns=["Latitudes", "Longitudes"])
df.to_csv(home_folder + 'ArcAtl_LatLon_H3_Res6.csv', index=False)
