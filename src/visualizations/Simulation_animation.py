import xarray as xr
import matplotlib.pyplot as plt
from datetime import timedelta
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.colors as color


home_folder = '/nethome/manra003/analysis/dispersion/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
mask_lon = mask_ds['glamf'].values
mask_lat = mask_ds['gphif'].values
mask_land = mask_ds['tmask'].values[:,0,:,:]

ds = xr.open_zarr(home_folder + 'simulations/Benguela_0pt0625_FULLm_Oct01-30_2015.zarr')
print(ds)

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['gray', 'whitesmoke'])

# region Animation
output_dt = timedelta(hours=6)
time_range = np.arange(np.nanmin(ds['time'].values),
                       np.nanmax(ds['time'].values) + np.timedelta64(output_dt),
                       output_dt)
print('Time_range: ', len(time_range))

# release locations
time_id = np.where(ds['time'] == time_range[0])
scatter = ax.scatter(ds['lon'].values[time_id], ds['lat'].values[time_id], s=1, c='sandybrown')

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(mask_lon[0, 999:2000, 1499:], mask_lat[0, 999:2000, 1499:], mask_land[0, 1000:2000, 1500:], cmap=colormap)

t = np.datetime_as_string(time_range[0], unit='m')
title = ax.set_title('Particles at z = 0 m and time = ' + t)


def animate(i):
    t = np.datetime_as_string(time_range[i], unit='m')

    time_id = np.where(ds['time'] == time_range[i])
    title.set_text('Particles at z = 0 m and time = ' + t)

    scatter.set_offsets(np.c_[ds['lon'].values[time_id], ds['lat'].values[time_id]])
    

size = len(time_range)
anim = FuncAnimation(fig, animate, frames=size, interval=200)
anim.save(home_folder + 'outputs/animations/Benguela_0pt0625_FULLm_Oct01-30_2015.mp4')
# endregion

print('animation saved')
