import xarray as xr
import numpy as np
import numpy as np
import dispersion_utils as utils
import os
import sys

home_folder = "/nethome/manra003/analysis/dispersion/simulations/"
output_folder = "/nethome/manra003/analysis/dispersion/outputs/"

min_lat, max_lat = -35, -15
min_lon, max_lon = 10, 20

months = [str(x).zfill(2) for x in range(1,13)]
years = np.arange(2010,2022)

dvm_types=['DVM', 'DVM_s','DVM_d']
compare_modes=['2D', '3D', 'sinking']

p_total = 104636    #total number of particles per simulation
t_days=100  #length of simulation in days

threshold_dist = 100    # km-  mesoscale

def get_separation_array(ds1, ds2, p_total, t_days):
    '''
    Method to get per particle pair distance between two types for each day opf the simulation
    parameters: ds1, ds2-two datasets; p_total-total number of particles; t_days-total number of outputs(per day) in the simulation 
    return: an array with horizontal separation distance between two types of particles for the same release.
    '''
    sep_array_1_2 = np.empty((p_total, t_days))
    sep_array_1_2[:] = np.nan
    for time_step in range(t_days):
        sep_array_1_2[:, time_step] = utils.distance(ds1['lon'][:, time_step], ds1['lat'][:, time_step], ds2['lon'][:, time_step], ds2['lat'][:, time_step])
    return sep_array_1_2


print("starting analysis")

for dvm_mode in dvm_types:
    for mode in compare_modes:
        distance_folder = output_folder + "Distances/{0}-{1}/".format(dvm_mode, mode)
        if not os.path.exists(distance_folder):
            os.makedirs(distance_folder)
        cdf_folder = output_folder + "CDF/{0}-{1}/".format(dvm_mode, mode)
        if not os.path.exists(cdf_folder):
            os.makedirs(cdf_folder)
        for year in years:
            all_sep_array = np.empty((p_total, t_days, len(months)))
            all_CDFs = np.zeros((len(months), t_days+1))
            print("Summary {0}-{1} for all months: {2}".format(dvm_mode, mode, year))

            # Get monthly CDF from 100 days simulation.
            for index, month in enumerate(months):
                ds1 = xr.open_zarr(home_folder + '{0}/{1}/Benguela_{0}_1ov32_641x_321yres_{1}-{2}_5z_100days.zarr'.format(dvm_mode, year, month))
                ds2 = xr.open_zarr(home_folder + '{0}/{1}/Benguela_{0}_1ov32_641x_321yres_{1}-{2}_5z_100days.zarr'.format(mode, year, month))

                # assert files are similar
                assert len(ds1.trajectory) == len(ds2.trajectory)
                assert len(ds1.obs) ==100 or len(ds1.obs) ==101
                assert len(ds2.obs) ==100 or len(ds2.obs) ==101  
                all_sep_array[:, :, index] = get_separation_array(ds1, ds2, p_total, t_days)
                all_CDFs[index, :], _ = utils.get_diff_CDF_PDF(all_sep_array[:, :, index], threshold_dist, t_days)

            np.save(distance_folder+"all_sep_array_{0}_{1}_{2}_tc{3}km.npy".format(dvm_mode, mode, year, threshold_dist), all_sep_array)
            np.save(cdf_folder+"all_CDF_{0}_{1}_{2}_tc{3}km.npy".format(dvm_mode, mode, year, threshold_dist), all_CDFs)
            sys.stdout.flush()