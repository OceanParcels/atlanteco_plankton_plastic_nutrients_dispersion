import pandas as pd
import h3
import numpy as np
import xarray as xr
from glob import glob
import adjacency_graph as ag
import connectivity_helper as ch

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
data_folder = '/Users/dmanral/Desktop/Analysis/UvA/'
hex_res = 3


def get_stationcode_hexes_mapping():
    stations = pd.read_csv(data_folder + 'Stations.csv', header=0, index_col=1)
    # sort stations in order of their latitudes
    sorted_stations = stations.sort_values(by=['Latitude'], ascending=False)
    return sorted_stations.index.values, sorted_stations['H3Id']


def main():
    domain_adjacency_file = 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz'

    temp_constraint_range = np.NaN
    min_accept_temp = np.NAN
    max_accept_temp = np.NAN
    prob_cutoff = np.NAN
    trahms2021 = False
    minimum_time_omalley2021 = False

    if minimum_time_omalley2021:
        t_ratio_file = home_folder + 'ProcessedTM/MinTRatio_FullAdjacency_csr.npz'
    else:
        t_ratio_file = None
    # Get all sorted stations- station codes between -100 and 20 Longitude
    stations_code, stations_hex = get_stationcode_hexes_mapping()

    master_hex_ids = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList.npy')
    # master_hex_ids = get_all_grids_hex_ids()
    # map station to master hex (some stations lie in the same hex- same connectivity)
    domain_mask = np.in1d(stations_hex, master_hex_ids)
    final_stations_code = stations_code[domain_mask]
    final_stations_hex = stations_hex[domain_mask]
    final_stations_master_indices = np.array([master_hex_ids.index(h) for h in final_stations_hex])
    # not use this as it sorts the results and return unique: np.in1d(master_hex_ids,stations_hex).nonzero()[0]

    # create graph
    if ~np.isnan(temp_constraint_range):
        print('Temp range: ', temp_constraint_range)
        atlantic_graph = ag.create_temp_range_graph(home_folder + domain_adjacency_file,
                                                    home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                    home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                    temp_constraint_range,
                                                    t_ratio_file)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        print('Min/Max Temp: ', min_accept_temp, max_accept_temp)
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + domain_adjacency_file,
                                                      home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                      home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                      min_accept_temp, max_accept_temp,
                                                      t_ratio_file)
    elif ~np.isnan(prob_cutoff):
        print('Probability Cutoff Value: ', prob_cutoff)
        atlantic_graph = ag.create_prob_filtered_graph(home_folder + domain_adjacency_file,
                                                       home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                       home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                       prob_cutoff)
    elif trahms2021:
        atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
                                              home_folder + 'ProcessedTM/Annual_SUM_MinTemperature_csr.npz',
                                              home_folder + 'ProcessedTM/Annual_SUM_MaxTemperature_csr.npz')

    else:
#         atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
#                                               home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
#                                               home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
#                                               t_ratio_file)
        atlantic_graph = ag.create_simple_graph(home_folder + domain_adjacency_file,
                                                t_ratio_file)

    # get min_T paths for all pairs- forward and backward- 2 d matrix.
    min_T_matrix = np.empty((len(final_stations_hex), len(final_stations_hex)))
    min_T_matrix[:] = np.NAN

    station_count = len(final_stations_master_indices)
    nnz_nan_count = 0
    zero_count = 0

    for i in range(station_count):
        s_idx = final_stations_master_indices[i]
        for j in range(i, station_count):
            d_idx = final_stations_master_indices[j]
            if s_idx != d_idx:
                f_path = ag.get_shortest_path(atlantic_graph, s_idx, d_idx)
                # f_path = ag.get_most_probable_path(atlantic_graph, s_idx, d_idx)
                # f_path = ag.minimum_time_path_malley2021(atlantic_graph, s_idx, d_idx)
                if f_path:
                    min_T_matrix[i][j] = len(f_path) - 1  # (time in months)
                    # min_T_matrix[i][j] = ag.get_time_from_most_probable_path(atlantic_graph, f_path)[-1]
                    nnz_nan_count += 1

                b_path = ag.get_shortest_path(atlantic_graph, d_idx, s_idx)
                # b_path = ag.get_most_probable_path(atlantic_graph, d_idx, s_idx)
                # b_path = ag.minimum_time_path_malley2021(atlantic_graph, d_idx, s_idx)
                if b_path:
                    min_T_matrix[j][i] = len(b_path) - 1  # (time in months)
                    # min_T_matrix[j][i] = ag.get_time_from_most_probable_path(atlantic_graph, b_path)[-1]
                    nnz_nan_count += 1
            else:
                if ag.check_if_edge_exists(atlantic_graph, s_idx, d_idx):
                    if i != j:
                        min_T_matrix[j][i] = 0
                        zero_count += 1
                    min_T_matrix[i][j] = 0
                    zero_count += 1

    print('maximum time: ', np.nanmax(min_T_matrix))
    print('Non zero-non NAN counter: ', nnz_nan_count)
    print('zero counter: ', zero_count)
    print("NAN count: ", np.count_nonzero(np.isnan(min_T_matrix)))
    print("non_zero count: ", np.count_nonzero(min_T_matrix))
    print("zero count: ", len(np.where(min_T_matrix == 0)[0]))
    assert nnz_nan_count + np.count_nonzero(np.isnan(min_T_matrix)) + len(np.where(min_T_matrix == 0)[0]) == \
           np.square(len(final_stations_hex))
    np.savez_compressed(
        data_folder + 'Stations_MinT_UvAconnectivity_nan.npz',
        codes=final_stations_code, matrix=min_T_matrix)


if __name__ == '__main__':
    main()