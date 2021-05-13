import xarray as xr
import matplotlib.pyplot as plt
from datetime import timedelta
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.colors as color
import pandas as pd

home_folder = '/Users/dmanral/Desktop/Analysis/AmazonOutflow/Task7_D/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

ds = xr.open_dataset(home_folder + 'Atlanteco_Romeo_July-Dec2018_OFF400km_z50m.nc')
print(ds)

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['gray', 'whitesmoke'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center

# Near Amazon river- 4-5 month simulations
ax.pcolormesh(x[0, 1599:2500, 99:1200], y[0, 1599:2500, 99:1200], c[0, 0, 1600:2500, 100:1200], cmap='Greens_r')
# plt.show()

# region Animation
output_dt = timedelta(days=5)
time_range = np.arange(np.nanmin(ds['time'].values),
                       np.nanmax(ds['time'].values) + np.timedelta64(output_dt),
                       output_dt)
print('Time_range: ', len(time_range))

# release locations
time_id = np.where(ds['time'] == time_range[0])
scatter = ax.scatter(ds['lon'].values[time_id], ds['lat'].values[time_id], s=1, c='sandybrown')

t = np.datetime_as_string(time_range[0], unit='m')
title = ax.set_title('Particles at z = 0 m and time = ' + t)


def animate(i):
    t = np.datetime_as_string(time_range[i], unit='m')

    time_id = np.where(ds['time'] == time_range[i])
    title.set_text('Particles at z = 50 m and time = ' + t)

    scatter.set_offsets(np.c_[ds['lon'].values[time_id], ds['lat'].values[time_id]])
    

size = len(time_range)
anim = FuncAnimation(fig, animate, frames=size, interval=200)
anim.save(home_folder + 'Romeo_July-Nov2018_OFF400km_z50m.mp4')
# endregion

print('animation saved')
