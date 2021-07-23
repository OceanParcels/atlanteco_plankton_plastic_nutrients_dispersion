import xarray as xr
import numpy as np
import pandas as pd
import h3
from scipy.sparse import coo_matrix

ds = xr.open_dataset('/Users/dmanral/Desktop/Analysis/TARA/Task4/tara_res5_01/FullTara_Res5_TS_May2014_dt600.nc')
hex_res = 3

hex_t0 = [h3.geo_to_h3(y, x, hex_res) for x, y in zip(ds['lon'][:, 0].values, ds['lat'][:, 0].values)]
hex_t1 = [h3.geo_to_h3(y, x, hex_res) for x, y in zip(ds['lon'][:, -1].values, ds['lat'][:, -1].values)]

hex_indices = np.unique(hex_t0)
n_total = len(hex_indices)
# mask hex ids that are new
hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, 'new')

# mask hex ids in hex_t1 that were deleted during the simulation
hex_t1_new = np.where(ds['time'][:, -1].values < np.max(ds['time'][:, -1].values), 'deleted', hex_t1_new)

hex_indices = np.append(hex_indices, ('new', 'deleted'))
mat_indices = np.arange(0, len(hex_indices))
map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

rows = map_h3_to_mat[hex_t0].values
cols = map_h3_to_mat[hex_t1_new].values
data = np.ones((len(hex_t0)))

t_matrix = coo_matrix((data, (rows, cols)), shape=(n_total, n_total + 2))
t_matrix.sum_duplicates()
