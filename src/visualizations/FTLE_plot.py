"""
source: https://github.com/LauraGomezNavarro/OceanParcels_Lyapunov/blob/main/code/FTLE_func_test.py
Author: Laura Gomez Navarro
"""
from math import sin, cos, sqrt, atan2, radians
import numpy as np
import xarray as xr
import numpy.linalg as LA
from datetime import timedelta, date
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib import colors
import sys

def ftle_brunton_2009(J, Td):  # http://cwrowley.princeton.edu/papers/BruntonChaos09.pdf
    D = np.dot(np.transpose(J), J)  # Cauchy–Green strain tensor
    lamda = LA.eigvals(D)
    lam_max = max(lamda)
    ftle = (1 / Td) * np.log(np.sqrt(lam_max))
    return ftle


def dist_pairs_km(inlon1, inlon2, inlat1, inlat2):
    """
    source: https://stackoverflow.com/questions/19412462/getting-distance-between-two-points-based-on-latitude-longitude

    """
    # approximate radius of earth in km
    R = 6373.0

    lon1 = radians(inlon1)
    lat1 = radians(inlat1)
    lon2 = radians(inlon2)
    lat2 = radians(inlat2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c

    return distance


args = sys.argv
assert len(args) == 4
year = np.int32(args[1])
month_num = np.int32(args[2])
mon_name = date(1900, month_num, 1).strftime('%b')
Td = np.int32(args[3])

home_folder = '/nethome/manra003/analysis/dispersion/'
output_folder = home_folder+ 'outputs/ftle/'

ds = xr.open_zarr(home_folder + 'simulations/Benguela_0625_401x257_{0}01-31_{1}.zarr'.format(mon_name, year))

output_dt = timedelta(days=1)  # from the simulation

out_index = Td # output everyday
time_range = np.arange(np.nanmin(ds['time'].values),
                       np.nanmax(ds['time'].values) + np.timedelta64(output_dt),
                       output_dt)
print(time_range[0], time_range[out_index])

coords = np.load(home_folder + 'Benguela_0625_401x257_release_points.npz')

# grid_lons, grid_lats = np.meshgrid(coords['Longitude'], coords['Latitude'])
grid_lons, grid_lats=ds['lon'][:, 0], ds['lat'][:, 0]

# initial position
x0 = np.reshape(ds['lon'][:, 0].data, (coords['Longitude'].shape[0], coords['Longitude'].shape[1]))
y0 = np.reshape(ds['lat'][:, 0].data, (coords['Longitude'].shape[0], coords['Longitude'].shape[1]))

# final position

x1 = np.reshape(ds['lon'][:, out_index].data, (coords['Longitude'].shape[0], coords['Longitude'].shape[1]))
y1 = np.reshape(ds['lat'][:, out_index].data, (coords['Longitude'].shape[0], coords['Longitude'].shape[1]))

H = x0.shape[0]
L = x1.shape[1]

FTLE_f = np.ones_like(np.asarray(x0))
FTLE_f[:,:] = np.NaN

J = np.empty([2, 2], float)

# 1, H-1 --> to ignore bordersx for now
for i in range(1, H - 1):  # 0, H-2
    for j in range(1, L - 1):  # 0, L-2
        J [:,:] = np.NaN
        J[0][0] = dist_pairs_km(x1[i, j], x1[i - 1, j], y1[i, j], y1[i - 1, j]) / dist_pairs_km(x0[i, j], x0[i - 1, j],
                                                                                                y0[i, j], y0[i - 1, j])
        J[0][1] = dist_pairs_km(x1[i, j], x1[i, j - 1], y1[i, j], y1[i, j - 1]) / dist_pairs_km(x0[i, j], x0[i, j - 1],
                                                                                                y0[i, j], y0[i, j - 1])
        J[1][0] = dist_pairs_km(x1[i, j], x1[i, j + 1], y1[i, j], y1[i, j + 1]) / dist_pairs_km(x0[i, j], x0[i, j + 1],
                                                                                                y0[i, j], y0[i, j + 1])
        J[1][1] = dist_pairs_km(x1[i, j], x1[i + 1, j], y1[i, j], y1[i + 1, j]) / dist_pairs_km(x0[i, j], x0[i + 1, j],
                                                                                                y0[i, j], y0[i + 1, j])

        if np.isnan(J).any():
            continue
        elif np.all((J == 1.0)):  # identify all land points- static points: No dispersion?
            continue
        else:
            f_value = ftle_brunton_2009(J, Td)
            FTLE_f[i][j] = f_value
print(np.nanmin(FTLE_f), np.nanmax(FTLE_f))

savename = output_folder + 'FTLE_BU_0625_401x257_{0}_{1}_{2}D.npz'.format(mon_name, year, Td)
np.savez(savename, FTLE_f=FTLE_f)

model_mask_file = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
mask_lon = mask_ds['glamf'].values
mask_lat = mask_ds['gphif'].values
mask_land = mask_ds['tmask'].values[:,0,:,:]

custom_size=10
fig = plt.figure(figsize=(12,8), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': custom_size, 'color': 'k'}
gl.ylabel_style = {'size': custom_size, 'color': 'k'}
ax.set_xlim(5, 21)
ax.set_ylim(-40, -10)
colormap = colors.ListedColormap(['gainsboro', 'white'])
ax.pcolormesh(mask_lon[0, 1249:1750, 1499:], mask_lat[0, 1249:1750, 1499:], mask_land[0, 1250:1750, 1500:], cmap=colormap)

plt.scatter(coords['Longitude'], coords['Latitude'], c=FTLE_f, cmap='RdBu_r', s=1)
plt.title('FLTE computation for {0} {1} after {4} days of release\n Minimum: {2} and Maximum: {3}'.format(mon_name, year, np.round(np.nanmin(FTLE_f),2), np.round(np.nanmax(FTLE_f),2), Td))
cbar = plt.colorbar()
cbar.set_label("FTLE (1/days)")
plt.clim(-0.4, 0.4)

plt.savefig(output_folder + "Benguela_0625_401x257_{0}_{1}_{2}D.jpeg".format(mon_name, year, Td))
