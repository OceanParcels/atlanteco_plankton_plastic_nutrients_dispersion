import numpy as np
from glob import glob
import xarray as xr
import h3
import pandas as pd
from time import time


def get_hex_id(lons, lats):
    return [h3.geo_to_h3(y, x, hex_resolution) for x, y in zip(lons, lats)]


def get_particle_info(ds, t):
    time_id = np.where(ds['time'] == t)
    particle_id = time_id[0]
    all_lats = ds['lat'].values[time_id]
    all_lons = ds['lon'].values[time_id]

    return particle_id, get_hex_id(all_lons, all_lats)


def verify_date(date, day, month, yr):
    d = pd.to_datetime(date)
    if d.day == day and d.month == month and d.year == yr:
        return date
    raise ValueError


home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
# seed_points = pd.read_csv(home_folder + 'Nemo_H3Release_LatLon_Res5_bs36.csv')
# no_particles = len(seed_points['Longitudes'])

# this doesnt seem to match with the order of elements from the simulation,
# therefore using second approach- read from a netcdf run
# master_hexId = seed_points['Res3_HexId'].unique()
temp_file = home_folder + '/tara_res5_01/FullTara_Res5_TS_Jun2014_dt600.nc'
temp_ds = xr.open_dataset(temp_file).load()
hex_resolution = 3
t_particleId, t_hexId = get_particle_info(temp_ds, np.min(temp_ds['time'][:, 0].values))

master_hexId = pd.unique(t_hexId)
no_particles = len(t_particleId)
master_particleId = np.arange(0, no_particles, 1)
assert (master_particleId == t_particleId).all()
temp_ds.close()

run_for_months = 12
total_no_hex = len(master_hexId)
# two extra columns added: new and deleted state
full_transition_matrix = np.zeros((run_for_months, total_no_hex, total_no_hex + 2))

months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
years = np.arange(2009, 2019, 1)
# currently 10 for each year simulation, will change with more simulations
per_month_simulations = 10

# region: transition_matrix
for mon_index, mon in zip(np.arange(0, run_for_months, 1), months[:run_for_months]):
    # files from all the years for that month
    monthly_files = sorted(glob(home_folder + '/tara_res5_01/FullTara_Res5_TS_{0}*_dt600.nc'.format(mon)))
    for file, year in zip(monthly_files, years):
        ds = xr.open_dataset(file).load()
        t_min, t_max = (np.min(ds['time'][:, 0].values), np.max(ds['time'][:, 1].values))
        ini_particleId, ini_hexId = get_particle_info(ds, t_min)
        final_particleId, final_hexId = get_particle_info(ds, t_max)
        # essential check to ensure that the same order of master list is used
        assert np.array_equal(pd.unique(ini_hexId), master_hexId)
        ds.close()
        st = time()
        print('timer started: {0} {1}'.format(mon, year))
        # create Transition Matrix for particles that were not deleted
        for i in range(len(final_particleId)):
            source_index = np.where(master_hexId == ini_hexId[final_particleId[i]])[0]
            assert len(source_index) == 1
            destination_index = np.where(master_hexId == final_hexId[i])[0]

            if len(destination_index) == 0:
                # update the value in the second-last column- particle moved to a new hex in model domain
                full_transition_matrix[mon_index, source_index, total_no_hex] += 1
            else:
                assert len(destination_index) == 1
                full_transition_matrix[mon_index, source_index, destination_index] += 1
        ed = time()
        print('TM- not deleted: ', (ed - st) / 60)
        # update deleted particles- last column
        deleted_particleId = np.setdiff1d(master_particleId, final_particleId)
        for d in deleted_particleId:
            deleted_index = np.where(master_hexId == ini_hexId[d])[0]
            assert len(deleted_index) == 1
            full_transition_matrix[mon_index, deleted_index, total_no_hex + 1] += 1
        ed2 = time()
        print('TM- deleted: ', (ed2 - ed) / 60)

    # month wise statss
    mon_sum = np.sum(full_transition_matrix[mon_index])
    assert mon_sum == no_particles * per_month_simulations
    print("-------------------------------\nMonth/Year: %s/%d- \n Sum: %f" % (mon, year, mon_sum))
    print("New count from matrix: ", np.sum(full_transition_matrix[mon_index, :, total_no_hex]))
    print("Deleted count from matrix: ", np.sum(full_transition_matrix[mon_index, :, total_no_hex + 1]))
    assert (np.sum(full_transition_matrix[mon_index, :, :], axis=2) == per_month_simulations).all()
    static = np.where(full_transition_matrix[mon_index, :, :] == per_month_simulations)
    # print(static[1])
    # print(static[2])
    print("Maximum Sum 10: ", len(static[0]))
    cnt_same = np.count_nonzero(
        [full_transition_matrix[mon_index, i, i] == per_month_simulations for i in range(total_no_hex)])
    print("Particles in same grid: ", len(cnt_same))
# endregion

# get the transition probability-matrix: sum of each row is 1
probability_matrix = full_transition_matrix / per_month_simulations
assert (np.round(np.sum(probability_matrix, axis=2), 4) == 1).all()
# assert (np.sum(probability_matrix, axis=2) == 1).all()# not working exactly= min 0.9999999..., max =1

# Matrix Analysis 3: Get average of connectivity over all the months.
avg_adjacency_matrix = np.average(probability_matrix, axis=0)
np.save(home_folder + 'Annual_Average_adjacency_matrix_Res5.npy', avg_adjacency_matrix)
print(max(avg_adjacency_matrix[0]))
print(max(avg_adjacency_matrix[1]))
