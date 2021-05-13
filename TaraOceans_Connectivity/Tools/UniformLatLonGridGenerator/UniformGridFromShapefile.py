# Shapefile source: https://www.marineregions.org/gazetteer.php?p=details&id=1902

import shapefile as shp
from shapely.ops import unary_union
from shapely.geometry import shape as shape_geom
from shapely.geometry import Point
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import time
import pandas as pd

home_folder = '/Tools/LatLonGridGenerator/'
file = home_folder + 'iho/iho.shp'

sf = shp.Reader(file)
print(sf)

poly_list = sf.shapes()

t1 = time.time()
# dissolve multiple polygons to get one outer bound
union = unary_union([shape_geom(s.__geo_interface__) for s in poly_list])  # returns shapely type polygon
print("Time to get the outer Polygon: ", time.time() - t1)

# Full atlantic
min_lat, min_lon = -60, -100
max_lat, max_lon = 70, 20

grid_size = 1
lats = np.arange(min_lat, max_lat + grid_size, grid_size)
lons = np.arange(min_lon, max_lon + grid_size, grid_size)

yy, xx = np.meshgrid(lons, lats)

fig = plt.figure()
ax = plt.axes()
ax.scatter(yy, xx, c='r', s=0.1)

x, y = union.exterior.xy

plt.scatter(x, y, s=0.1)

print("number of lons: ", len(lons))
print("number of lats: ", len(lats))

t1 = time.time()
Points = [Point(x, y) for x, y in product(lons, lats)]
print("Total set of points: %d for Grid size of: %f", len(Points), grid_size)
Points = unary_union(Points)

# intersection over within or contains- this one return points(MultiPoint)  that lie within the polygon
result = Points.intersection(union)
print("Time to get the points within the outer polygon: ", time.time() - t1)

print("Points within the Atlantic Ocean: ", len(result))

final_lats = [a.y for a in result]
final_lons = [a.x for a in result]

ax.scatter(final_lons, final_lats, s=3, c='g')
plt.show()

coords = {"Latitudes": final_lats, "Longitudes": final_lons}
df = pd.DataFrame(coords, columns=["Latitudes", "Longitudes"])
df.to_csv(r'/scratch/dmanral/Atlantic_LatLons_1degree_mesh.csv', index=False)
