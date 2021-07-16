import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as color
import numpy as np
import pandas as pd
import h3

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

ds = xr.open_dataset(home_folder + 'Test_FullTara_Res6_Jan2018_landpoints.nc')
print(ds)
total_lats = ds['lat'][:, 0]
total_lons = ds['lon'][:, 0]

static_pts = ds.where((ds['lat'][:, 0] == ds['lat'][:, -1]) & (ds['lon'][:, 0] == ds['lon'][:, -1]), drop=True)
s_lats = static_pts['lat'][:, 0]
s_lons = static_pts['lon'][:, 0]
print(len(s_lons))

# remove these points from the seed points list
seed_lats_bool = np.in1d(total_lats, s_lats)
seed_lons_bool = np.in1d(total_lons, s_lons)
assert np.array_equal(seed_lats_bool, seed_lons_bool)

index = np.where(seed_lons_bool == False)
seed_lats = np.array(total_lats.loc[index[0]])
seed_lons = np.array(total_lons.loc[index[0]])
print(len(seed_lons))

# Get Points that were deleted from the simulation at time=0, as they were out of the model domain
# Sample format listed below
extra_points = [[[-78.374268, -70.189794],
                [-78.376736, -67.233323],
                [-78.379466, -64.001345],
                [-78.376418, -63.462016]]

extra_lats = [a[0] for a in extra_points]
extra_lons = [a[1] for a in extra_points]

# remove these points from the seed points list
seed_lats_bool = np.in1d(np.round(seed_lats, 6), extra_lats)
seed_lons_bool = np.in1d(np.round(seed_lons, 6), extra_lons)
# assert np.array_equal(seed_lats_bool, seed_lons_bool)

index = np.where(seed_lons_bool == False)
# seed_lats = np.array(seed_lats['Latitudes'].loc[index])
# seed_lons = np.array(seed_lons['Longitudes'].loc[index])

seed_lats = np.array(seed_lats[index])
seed_lons = np.array(seed_lons[index])

print(len(seed_lons))

# Add locations for Stations outside the hexagons- Do not use the centroids for these cells
# OA247
# station_lats = np.array([51.55385])
# station_lons = np.array([-8.373833333])

# full_lats = np.append(seed_lats, station_lats)
# full_lons = np.append(seed_lons, station_lons)
# print(full_lats.shape, full_lons.shape)

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['gainsboro', 'whitesmoke'])

ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
# Plot points to delete
# ax.scatter(s_lons, s_lats, s=1, c='r')

ax.scatter(seed_lons, seed_lats, s=0.5, c='g')
# ax.scatter(station_lons, station_lats, s=0.5, c='r')

plt.show()

hex_resolution = 3
res3_hex_list = [h3.geo_to_h3(y, x, hex_resolution) for x, y in zip(seed_lons, seed_lats)]

unique_hex = np.unique(res3_hex_list)
print("total number of unique release hex grids", unique_hex.shape)

coords = {"Latitudes": seed_lats, "Longitudes": seed_lons, "Res3_HexId": res3_hex_list}
df = pd.DataFrame(coords, columns=["Latitudes", "Longitudes", "Res3_HexId"])
df.to_csv(home_folder + 'Nemo_H3Release_LatLon_Res6.csv', index=False)
