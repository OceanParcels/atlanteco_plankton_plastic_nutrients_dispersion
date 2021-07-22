import numpy as np
from glob import glob
import xarray as xr
import h3
import pandas as pd


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


home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task3_TM/'
master_particleId = np.arange(0, 8114, 1)
seed_points = pd.read_csv(home_folder + 'Nemo_H3Release_LatLon.csv')
master_hexId = seed_points['Res3_HexId'].unique()

run_for_months = 12
total_no_particles = len(master_hexId)
# two extra columns added: new and deleted state
full_transition_matrix = np.zeros((run_for_months, total_no_particles, total_no_particles + 2))

months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
years = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]
# currently 10 for each year simulation, will change with more simulations
per_month_simulations = 10
hex_resolution = 3

# region: transition_matrix
for mon in months[:run_for_months]:
    # files from all the years for that month
    monthly_files = sorted(glob(home_folder + '/tara_data/FullAtlantic_2D_01{0}*_1month.nc'.format(mon)))
    mon_index = np.where(months == mon)[0]

    for file, year in zip(monthly_files, years):
        ds = xr.open_dataset(file).load()

        # ensure start date and end date are as expected
        t_min = verify_date(np.nanmin(ds['time'].values), 1, mon_index + 1, year)
        if mon_index < 11:
            t_max = verify_date(np.nanmax(ds['time'].values), 1, mon_index + 2, year)
        else:
            if year != 2018:
                t_max = verify_date(np.nanmax(ds['time'].values), 1, 1, year + 1)
            else:
                t_max = verify_date(np.nanmax(ds['time'].values), 31, mon_index + 1, year)

        final_particleId, final_hexId = get_particle_info(ds, t_max)

        # verify that the order of the initial hex_ids from the file is same as in the master hex_ids
        ini_particleId, ini_hexId = get_particle_info(ds, t_min)
        assert np.array_equal(ini_hexId, master_hexId)

        ds.close()

        # create Transition Matrix for particles that were not deleted
        for i in range(len(final_particleId)):
            f_index = np.where(master_hexId == final_hexId[i])[0]
            if len(f_index) == 0:
                # update the value in the last column
                full_transition_matrix[mon_index, final_particleId[i], total_no_particles] += 1
            else:
                full_transition_matrix[mon_index, final_particleId[i], f_index] += 1

        # update deleted particles count- second-last column
        deleted = np.setdiff1d(master_particleId, final_particleId)
        full_transition_matrix[mon_index, deleted, total_no_particles + 1] += 1

    # month wise stats
    mon_sum = np.sum(full_transition_matrix[mon_index])
    assert mon_sum == total_no_particles * per_month_simulations
    print("-------------------------------\nMonth/Year: %s/%d- \n Sum: %f" % (mon, year, mon_sum))
    print("New count from matrix: ", np.sum(full_transition_matrix[mon_index, :, total_no_particles]))
    print("Deleted count from matrix: ", np.sum(full_transition_matrix[mon_index, :, total_no_particles + 1]))
    assert (np.sum(full_transition_matrix[mon_index, :, :], axis=2) == per_month_simulations).all()
    static = np.where(full_transition_matrix[mon_index, :, :] == 10)
    print(static[1])
    print(static[2])
    print("Maximum Sum 10: ", len(static[0]))
    cnt_same = np.count_nonzero([full_transition_matrix[mon_index, i, i] == 10 for i in range(total_no_particles)])
    print("Particles in same grid: ", cnt_same)
# endregion

# get the transition probability-matrix: sum of each row is 1
probability_matrix = full_transition_matrix / per_month_simulations
assert (np.round(np.sum(probability_matrix, axis=2), 4) == 1).all()
# assert (np.sum(probability_matrix, axis=2) == 1).all()# not working exactly= min 0.9999999..., max =1

stations = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1)
atlantic_lon_index = np.where(np.logical_and(stations['Longitude'] >= -100, stations['Longitude'] <= 20))

atlantic_lon = np.take(stations['Longitude'], atlantic_lon_index[0])
atlantic_lat = np.take(stations['Latitude'], atlantic_lon_index[0])

stations_hex = get_hex_id(atlantic_lon, atlantic_lat)

atlantic_hex, ind_1, ind2 = np.intersect1d(stations_hex, master_hexId, return_indices=True)
# within atlantic(104 unique hex), only 3 hex Ids have 2 stations mapped to each: /task1/Stations_in_same_H3Grids.png

# Matrix Analysis 3: Get average of connectivity over all the months.
avg_adjacency_matrix = np.average(probability_matrix, axis=0)
np.save(home_folder + 'avg_adjacency_matrix.npy', avg_adjacency_matrix)

# temp = avg_adjacency_matrix[ind2, :]
# avg_adMatrix_stations = temp[:, ind2]

# endregion

print(max(avg_adjacency_matrix[0]))
print(max(avg_adjacency_matrix[1]))
# print(advection_matrix)