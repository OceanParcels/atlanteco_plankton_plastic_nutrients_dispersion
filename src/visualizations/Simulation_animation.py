import xarray as xr
import matplotlib.pyplot as plt
from datetime import timedelta, date
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.colors as color
import sys
import cartopy.crs as ccrs

args = sys.argv
assert len(args) == 5
year = np.int32(args[1])
month_num = np.int32(args[2])
exact_depth = args[3]
release_depth = int(np.float32(exact_depth))
rk_mode = args[4]

home_folder = '/nethome/manra003/analysis/dispersion/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
mask_lon = mask_ds['glamf'].values
mask_lat = mask_ds['gphif'].values
mask_land = mask_ds['tmask'].values[:, 0, :, :]
file = home_folder + "simulations/{0}/{1}/Benguela_{0}_1ov32_641x_321yres_{1}-{2}_{3}z_{4}days.zarr".format(
    rk_mode, year, str(month_num).zfill(2), release_depth, 100)
print(file)
ds = xr.open_zarr(file)

print(ds)

custom_size = 10

fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': custom_size, 'color': 'k'}
gl.ylabel_style = {'size': custom_size, 'color': 'k'}

colormap = color.ListedColormap(['gray', 'whitesmoke'])

# region Animation
output_dt = timedelta(days=1)
time_range = np.arange(np.nanmin(ds['time'].values),
                       np.nanmax(ds['time'].values) +
                       np.timedelta64(output_dt),
                       output_dt)
print('Time_range: ', len(time_range))

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(mask_lon[0, :2000, 605:], mask_lat[0, :2000,
              605:], mask_land[0, 1:2000, 606:], cmap=colormap)

# release locations
time_id = np.where(ds['time'] == time_range[0])
scatter = ax.scatter(ds['lon'].values[time_id],
                     ds['lat'].values[time_id], s=1, c='sandybrown')

t = np.datetime_as_string(time_range[0], unit='m')
title = ax.set_title('Particles at z = ' + exact_depth + 'm and time = ' + t)
ax.set_xlim(3, 21)
ax.set_ylim(-42, -13)


def animate(i):
    t = np.datetime_as_string(time_range[i], unit='m')

    time_id = np.where(ds['time'] == time_range[i])
    title.set_text('Particles at z = ' + exact_depth + 'm and time = ' + t)

    scatter.set_offsets(
        np.c_[ds['lon'].values[time_id], ds['lat'].values[time_id]])


size = len(time_range)
anim = FuncAnimation(fig, animate, frames=size, interval=200)
anim.save(home_folder + 'outputs/animations/Benguela_{0}_1ov32_641x_321yres_{1}-{2}_{3}z_{4}days.mp4'.format(
    rk_mode, year, str(month_num).zfill(2), release_depth, 100))

# endregion

print('animation saved')
