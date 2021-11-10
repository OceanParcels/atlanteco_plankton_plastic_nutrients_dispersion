import numpy as np
import pandas as pd
import xarray as xr
import h3
import AdjacencyGraph as ag
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

hex_res = 3
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'


def plot_paths(x, y, c, master_hex_ids, forward_path, backward_path, f_time_laps, b_time_laps, s_code, d_code,
               path_to_compute):
    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    ax.set_title('{0} Path and Time estimate between station {1} and {2}'.format(path_to_compute, s_code, d_code))

    # https://stackoverflow.com/questions/48520393/filling-shapefile-polygons-with-a-color-in-matplotlib
    def plot_path(time_laps, path, cmap, s, d):
        time_laps = np.append(0, time_laps)
        norm = plt.Normalize(vmin=0, vmax=time_laps[-1])
        patches = []

        for i in range(len(path)):
            color = cmap(norm(time_laps[i]))
            polygons = h3.h3_set_to_multi_polygon([master_hex_ids[path[i]]], geo_json=True)
            patches.append(Polygon(polygons[0][0], True, color=color))

        pc = PatchCollection(patches, match_original=True, edgecolor=None, linewidths=None, zorder=1)
        ax.add_collection(pc)
        sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm1.set_array(time_laps)
        cb1 = fig.colorbar(sm1, ax=ax)
        cb1.set_label('{0} to {1} path time laps in years'.format(s, d))

    plot_path(f_time_laps, forward_path, plt.cm.winter, s_code, d_code)
    plot_path(b_time_laps, backward_path, plt.cm.Wistia.reversed(), d_code, s_code)
    plt.show()


def plot_shortest_paths_subset(x, y, c, master_hex_ids, f_paths, b_paths, s_code, d_code):
    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    ax.set_title('Sample of Shortest Possible Paths between station {0} and {1}'.format(s_code, d_code))

    def plot_path(path, color):
        centers = [h3.h3_to_geo(master_hex_ids[ind]) for ind in path]
        lats = [x1[0] for x1 in centers]
        lons = [x1[1] for x1 in centers]
        ax.plot(lons, lats, color=color, linestyle='solid', alpha=0.3)
        ax.scatter(lons, lats, color=color, s=0.4)

    [plot_path(f_paths[i], 'blue') for i in range(len(f_paths))]
    [plot_path(b_paths[i], 'red') for i in range(len(b_paths))]
    plt.show()


def load_mask_file():
    mask_ds = xr.open_dataset(home_folder + 'GLOB16L98_mesh_mask_atlantic.nc', decode_times=False).load()
    x = mask_ds['glamf']
    y = mask_ds['gphif']
    c = mask_ds['tmask'][:]
    return x, y, c


def get_station_hexes_from_code(s_code, d_code):
    stations = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1, index_col=0)

    def get_station_hexid(code):
        try:
            lat, lon = stations.loc[code]['Latitude'], stations.loc[code]['Longitude']
            return h3.geo_to_h3(lat, lon, hex_res)
        except KeyError:
            print("Incorrect source/destination code provided. Recheck values")
            raise

    return get_station_hexid(s_code), get_station_hexid(d_code)


def get_all_grids_hex_ids():
    # hexId and matrix index mapping from a random simulation output file from the output dataset
    with xr.open_dataset(home_folder + '/tara_res5_01/FullTara_Res5_TS_Jun2014_dt600.nc').load() as temp_ds:
        master_all_hex_t0 = np.array(
            [h3.geo_to_h3(y, x, hex_res) for x, y in zip(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values)])
        master_uni_hex = np.unique(master_all_hex_t0)
        return master_uni_hex.tolist()


def main():
    source_code = '66SUR'
    destination_code = 'OA002'
    temp_constraint_range = np.NAN
    min_accept_temp = -1
    max_accept_temp = 3

    s_hex, d_hex = get_station_hexes_from_code(source_code, destination_code)
    master_grids_list = get_all_grids_hex_ids()
    try:
        s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
    except KeyError:
        print("Source/destination code not present in the domain. Recheck values")
        raise

    # graph where the min-max 'temperature range' b/w grids is restricted
    if ~np.isnan(temp_constraint_range):
        atlantic_graph = ag.create_temp_range_graph(home_folder + 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz',
                                                    home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                    home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                    temp_constraint_range)
    # graph where connections with min and max temperature values are within the min-max values provided by the user
    # (eg. a thermal range for a species)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz',
                                                      home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                      home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                      min_accept_temp, max_accept_temp)
    # simple graph, without any temperature boundaries
    else:
        atlantic_graph = ag.create_simple_graph(home_folder + 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz')
    print('Graph ready')

    # region: get paths and times
    mask_lons, mask_lats, mask_value = load_mask_file()
    path_to_compute = 'Most Likely'
    forward_path = ag.get_most_probable_path(atlantic_graph, s_index, d_index)
    f_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, forward_path)

    backward_path = ag.get_most_probable_path(atlantic_graph, d_index, s_index)
    b_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, backward_path)
    plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
               b_time_laps, source_code, destination_code, path_to_compute)

    path_to_compute = 'Shortest (Mininum Time)'
    forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
    f_time_laps = np.arange(1, len(forward_path), 1) / 12
    backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
    b_time_laps = np.arange(1, len(backward_path), 1) / 12
    plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
               b_time_laps, source_code, destination_code, path_to_compute)

    path_to_compute = '100 Shortest Paths'
    forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index, len(forward_path))
    backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_index, s_index, len(backward_path))
    plot_shortest_paths_subset(mask_lons, mask_lats, mask_value, master_grids_list, forward_paths, backward_paths,
                               source_code, destination_code)
    # endregion


if __name__ == '__main__':
    main()