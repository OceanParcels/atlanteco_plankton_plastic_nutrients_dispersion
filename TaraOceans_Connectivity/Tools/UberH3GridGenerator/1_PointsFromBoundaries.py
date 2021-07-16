import matplotlib.pyplot as plt
import h3
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.colors as color

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'

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
colormap = color.ListedColormap(['gainsboro', 'whitesmoke'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

seed_point_resolution = 6

# For test region: Gulf of Mexico
# geoJSON = {'type': 'Polygon',
#            'coordinates': [[[0, -98], [0, -70], [30, -70], [30, -98]]]}
# For whole study region: Atlantic + Arctic + SO
geoJSON = {'type': 'Polygon',
           'coordinates': [[[-78.38, -100], [-78.38, 20], [88.2, 20], [88.2, -100]]]}
[[[-85, -100], [-85, 20], [90, 20], [90, -100]]]
# small patch of pacific
# geoJSON = {'type': 'Polygon',
#            'coordinates': [[[-60, -75], [-60, -67.2], [-56, -67.2], [-56, -75]]]}

hexagons = list(h3.polyfill(geoJSON, seed_point_resolution))
print("number of hexagons in the region with grid resolution 6:", len(hexagons))

hex_count = len(hexagons)
cen_lats = np.zeros(hex_count)
cen_lons = np.zeros(hex_count)

centroids = [h3.h3_to_geo(hex) for hex in hexagons]
print(len(centroids))

cen_lats = [c[0] for c in centroids]
cen_lons = [c[1] for c in centroids]

ax.scatter(cen_lons, cen_lats, s=1, c='r')

# # plot the hex for another resolution
# hex_list = [h3.geo_to_h3(y, x, hex_resolution) for x, y in zip(cen_lons, cen_lats)]
#
# unique_hex = np.unique(hex_list)
# print(unique_hex.shape)


# so_points = pd.read_csv(home_folder + 'SO_LatLon_H3_Res4.csv')
# ax.scatter(so_points['Longitudes'],so_points['Latitudes'], s=0.5, c='g')
#
# at_points = pd.read_csv(home_folder + 'Atlantic_Arctic_LatLon_H3_Res4.csv')
# ax.scatter(at_points['Longitudes'],at_points['Latitudes'], s=0.5, c='b')

plt.show()

coords = {"Latitudes": cen_lats, "Longitudes": cen_lons}
df = pd.DataFrame(coords, columns=["Latitudes", "Longitudes"])
df.to_csv(home_folder + 'ROMEO_Res6_Centroids.csv', index=False)
