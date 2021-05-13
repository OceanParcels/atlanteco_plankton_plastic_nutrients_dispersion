import xarray as xr
import matplotlib.pyplot as plt
from datetime import timedelta
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
import matplotlib.colors

fig = plt.figure()
ax = plt.axes()
home_folder = '/Users/dmanral/Desktop/Analysis/AmazonOutflow/Task7_D/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

ds1 = xr.open_dataset(home_folder + 'Atlanteco_Romeo_July-Dec2018_50km.nc')
ds2 = xr.open_dataset(home_folder + 'Atlanteco_Romeo_July-Dec2018_100km.nc')
ds3 = xr.open_dataset(home_folder + 'Atlanteco_Romeo_July-Nov2018_OFF250km_z0m.nc')
ds4 = xr.open_dataset(home_folder + 'Atlanteco_Romeo_July-Nov2018_OFF400km_z0m.nc')

# remove the first row and first column from the glamf/gphif to access points enclosed in the center

colormap = matplotlib.colors.ListedColormap(['gray', 'whitesmoke'])

# Near Amazon river- 4-5 month simulations
ax.pcolormesh(x[0, 1699:2300, 99:1000], y[0, 1699:2300, 99:1000], c[0, 0, 1700:2300, 100:1000], cmap=colormap)


# region Animation
output_dt = timedelta(days=1)
time_range = np.arange(np.nanmin(ds3['time'].values),
                       np.nanmax(ds3['time'].values) + np.timedelta64(output_dt),
                       output_dt)
print('Time_range: ', len(time_range))

time_id = np.where(ds3['time'] == time_range[0])
custom_lines = [Line2D([0], [0], color='tomato', lw=4),
                Line2D([0], [0], color='gold', lw=4),
                Line2D([0], [0], color='limegreen', lw=4),
                Line2D([0], [0], color='royalblue', lw=4)]

scatter1 = ax.scatter(ds1['lon'].values[time_id], ds1['lat'].values[time_id], s=0.2, marker="o", c='tomato')
scatter2 = ax.scatter(ds2['lon'].values[time_id], ds2['lat'].values[time_id], s=0.2, alpha=0.7, marker="o", c='gold')
scatter3 = ax.scatter(ds3['lon'].values[time_id], ds3['lat'].values[time_id], s=0.2, alpha=0.7, marker="o",
                      c='limegreen')
scatter4 = ax.scatter(ds3['lon'].values[time_id], ds3['lat'].values[time_id], s=0.2, alpha=0.7, marker="o",
                      c='royalblue')

ax.legend(custom_lines, ['50 km', '100 km', '250 km', '400km'])

t = np.datetime_as_string(time_range[0], unit='m')
title = ax.set_title('Particles at z = 0 m  and time = ' + t)


def animate(i):
    t = np.datetime_as_string(time_range[i], unit='m')
    title.set_text('Particles at z = 0 m  and time = ' + t)

    time_id = np.where(ds3['time'] == time_range[i])
    scatter1.set_offsets(np.c_[ds1['lon'].values[time_id], ds1['lat'].values[time_id]])
    scatter2.set_offsets(np.c_[ds2['lon'].values[time_id], ds2['lat'].values[time_id]])
    scatter3.set_offsets(np.c_[ds3['lon'].values[time_id], ds3['lat'].values[time_id]])
    scatter4.set_offsets(np.c_[ds4['lon'].values[time_id], ds4['lat'].values[time_id]])


size = len(time_range)
anim = FuncAnimation(fig, animate, frames=size, interval=200)
anim.save(
    '/Users/dmanral/Desktop/Analysis/AmazonOutflow/Task7_D/Romeo_July-Nov2018_comparison_z0m.mp4')
# endregion

print('animation saved')
