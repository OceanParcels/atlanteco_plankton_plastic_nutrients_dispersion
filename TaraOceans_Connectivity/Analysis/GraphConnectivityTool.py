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


def plot_mlp_paths(master_hex_ids, forward_path, backward_path, f_time_laps, b_time_laps, s_code, d_code,
                   path_to_compute):
    # plot ocean model mask
    mask_ds = xr.open_dataset(home_folder + 'GLOB16L98_mesh_mask_atlantic.nc', decode_times=False).load()

    x = mask_ds['glamf']
    y = mask_ds['gphif']
    c = mask_ds['tmask'][:]

    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    ax.set_title(
        '{0} Path and Time estimate between station {1} and {2}'.format(path_to_compute, s_code, d_code))

    # https://stackoverflow.com/questions/48520393/filling-shapefile-polygons-with-a-color-in-matplotlib
    def plot_path(time_laps, path, cmap):
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
        cb1.set_label('{0} to {1} path time laps in years'.format(s_code, d_code))

    plot_path(f_time_laps, forward_path, plt.cm.winter)
    plot_path(b_time_laps, backward_path, plt.cm.Wistia.reversed())

    plt.show()


def plot_shortest_paths_subset(master_hex_ids, paths, s_code, d_code):
    # plot ocean model mask
    mask_ds = xr.open_dataset(home_folder + 'GLOB16L98_mesh_mask_atlantic.nc', decode_times=False).load()

    x = mask_ds['glamf']
    y = mask_ds['gphif']
    c = mask_ds['tmask'][:]

    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    ax.set_title(
        'Sample of Shortest Possible Paths between station {0} and {1}'.format(s_code, d_code))

    def plot_path(path, color):
        centers = [h3.h3_to_geo(master_hex_ids[ind]) for ind in path]
        lats = [x1[0] for x1 in centers]
        lons = [x1[1] for x1 in centers]
        ax.plot(lons, lats, color=color, linestyle='solid', alpha=0.3)
        ax.scatter(lons, lats, color=color, s=0.2)

    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)

    colors = get_cmap(len(paths))
    [plot_path(paths[i], colors(i)) for i in range(len(paths))]
    plt.show()


def main():
    # parameters from arguments?? pairs(66SUR,OA002)(84SUR,210SUR)
    source_code = '66SUR'
    destination_code = 'OA002'
    temp_constraint_range = 1.5
    path_to_compute = 'Shortest'

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

    s_hex, d_hex = get_station_hexes_from_code(source_code, destination_code)
    master_grids_list = get_all_grids_hex_ids()
    try:
        s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
    except KeyError:
        print("Source/destination code not present in the domain. Recheck values")
        raise

    if temp_constraint_range != 0:
        atlantic_graph = ag.create_temp_sal_graph(home_folder + 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz',
                                                  home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                  home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                  temp_constraint_range,
                                                  home_folder + 'ProcessedTM/Annual_Avg_MinSalinity_csr.npz',
                                                  home_folder + 'ProcessedTM/Annual_Avg_MaxSalinity_csr.npz')
    else:
        atlantic_graph = ag.create_simple_graph(home_folder + 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz')

    # region: get paths and times based on option selected (args)
    if path_to_compute == 'Most_Likely':
        forward_path = ag.get_most_probable_path(atlantic_graph, s_index, d_index)
        f_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, forward_path)

        backward_path = ag.get_most_probable_path(atlantic_graph, d_index, s_index)
        b_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, backward_path)
        plot_mlp_paths(master_grids_list, forward_path, backward_path, f_time_laps, b_time_laps, source_code,
                       destination_code, path_to_compute)

    elif path_to_compute == 'Shortest':
        forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
        f_time_laps = np.arange(1, len(forward_path), 1) / 12
        backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
        b_time_laps = np.arange(1, len(backward_path), 1) / 12
        plot_mlp_paths(master_grids_list, forward_path, backward_path, f_time_laps, b_time_laps, source_code,
                       destination_code, path_to_compute)

    elif path_to_compute == 'Shortest_Paths':
        forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
        paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index, len(forward_path))
        plot_shortest_paths_subset(master_grids_list, paths, source_code, destination_code)
        # backward_path = get_shortest_paths_subset(atlantic_graph, d_index, s_index)

    else:
        raise KeyError('correct method to compute paths not provided')
    # endregion


if __name__ == '__main__':
    main()
