import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.colors as clr
import h3
from shapely.geometry.polygon import Polygon

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'

files1 = pd.read_csv(home_folder + 'SO_LatLon_H3_Res6.csv')
print(files1.shape)
files2 = pd.read_csv(home_folder + 'ArcAtl_LatLon_H3_Res6.csv')
print(files2.shape)
files3 = pd.read_csv(home_folder + 'ROMEO_Res6_PacificCentroids.csv')
print(files3.shape)

lats = np.array(list(files1['Latitudes']) + list(files2['Latitudes']) + list(
    files3['Latitudes']))
lons = np.array(list(files1['Longitudes']) + list(files2['Longitudes']) + list(
    files3['Longitudes']))
print(len(lats), len(lons))
print(len(np.unique(lats)), len(np.unique(lons)))

# remove duplicates (8)
lats_index = np.unique(lats, return_index=True)[1]
lons_index = np.unique(lons, return_index=True)[1]
#
# assert np.equal(np.sort(lats_index), np.sort(lons_index)).all()
#
indices = np.sort(lons_index)
full_lats = lats[indices]
full_lons = lons[indices]
print(len(full_lats))

# # region: Get stations
# stations = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1)
#
# atlantic_lon_index = np.where(np.logical_and(stations['Longitude'] >= -100, stations['Longitude'] <= 20))
#
# atlantic_lon = np.take(stations['Longitude'], atlantic_lon_index[0])
# atlantic_lat = np.take(stations['Latitude'], atlantic_lon_index[0])
# # endregion
# region: Get the plot for model mask
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
# endregion

hex_resolution = 6
res5_hex_list = [h3.geo_to_h3(y, x, hex_resolution) for x, y in zip(full_lons, full_lats)]

unique_hex = np.unique(res5_hex_list)
print("unique release points", unique_hex.shape)

hex_resolution = 3
res3_hex_list = [h3.geo_to_h3(y, x, hex_resolution) for x, y in zip(full_lons, full_lats)]

unique_hex = np.unique(res3_hex_list)
print("total number of unique release hex grids", unique_hex.shape)

# def plot_hex(hex, c):
#     polygons = h3.h3_set_to_multi_polygon([hex], geo_json=False)
#     p = Polygon(polygons[0][0])
#     y, x = p.exterior.xy
#     ax.plot(x, y, color=c, linewidth=0.5)
#
#
# for h in unique_hex:
#     plot_hex(str(h), 'royalblue')

ax.scatter(full_lons, full_lats, s=0.2, c='g')

# ax.scatter(atlantic_lon, atlantic_lat, s=5, c='r')
plt.show()

coords = {"Latitudes": full_lats, "Longitudes": full_lons}
df = pd.DataFrame(coords, columns=["Latitudes", "Longitudes"])
df.to_csv(home_folder + 'Nemo_H3Release_LatLon_Res6_withLandPoints.csv', index=False)
