import xarray as xr
import matplotlib.pyplot as plt
from datetime import timedelta
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.colors as color
import pandas as pd

home_folder = '/Users/dmanral/Desktop/Analysis/AmazonOutflow/Task8_TM/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

ds = xr.open_dataset(home_folder + 'FullAtlantic_1degree_mesh_2months.nc')
print(ds)
stations = pd.read_excel(home_folder + 'AllStations.xls', header=1)
atlantic_lon_index = np.where(np.logical_and(stations['Longitude'] >= -100, stations['Longitude'] <= 20))

atlantic_lon = np.take(stations['Longitude'], atlantic_lon_index[0])
atlantic_lat = np.take(stations['Latitude'], atlantic_lon_index[0])

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['gray', 'whitesmoke'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center

# Near Amazon river- 4-5 month simulations
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

# region Animation
output_dt = timedelta(days=1)
time_range = np.arange(np.nanmin(ds['time'].values),
                       np.nanmax(ds['time'].values) + np.timedelta64(output_dt),
                       output_dt)
print('Time_range: ', len(time_range))

time_id = np.where(ds['time'] == time_range[0])
lons = ds['lon'].values[time_id]
scatter = ax.scatter(lons, ds['lat'].values[time_id], s=1, c='sandybrown')
scatter_stations = ax.scatter(atlantic_lon, atlantic_lat, s=1, c='r')

t = np.datetime_as_string(time_range[0], unit='m')
title = ax.set_title(str(len(lons)) + ' Particles at z = 0 m and time = ' + t)


def animate(i):
    t = np.datetime_as_string(time_range[i], unit='m')

    time_id = np.where(ds['time'] == time_range[i])
    lons = ds['lon'].values[time_id]
    title.set_text(str(len(lons)) + ' Particles at z = 0 m and time = ' + t)

    scatter.set_offsets(np.c_[lons, ds['lat'].values[time_id]])
    scatter_stations.set_offsets(np.c_[atlantic_lon, atlantic_lat])


size = len(time_range)
anim = FuncAnimation(fig, animate, frames=size, interval=200)
anim.save(home_folder + 'Full_atlantic_1deg_withStations.mp4')
# endregion

print('animation saved')
