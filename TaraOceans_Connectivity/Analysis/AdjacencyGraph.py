import graph_tool.all as gt
from scipy.sparse import load_npz, save_npz, csr_matrix
import numpy as np


def get_base_graph(v_count, nnz_index, weights):
    g = gt.Graph()
    g.add_vertex(v_count)
    g.add_edge_list(np.transpose(nnz_index))
    ew = g.new_edge_property('double')
    ew.a = -np.log(weights)
    g.ep['weight'] = ew

    eprob = g.new_edge_property('double')
    eprob.a = weights
    g.ep['probability'] = eprob

    return g


def add_temperature_properties(g, min_temp, max_temp):
    min_t = g.new_edge_property('double')
    min_t.a = min_temp
    g.ep['min_t'] = min_t
    max_t = g.new_edge_property('double')
    max_t.a = max_temp
    g.ep['max_t'] = max_t
    return g


def add_minimum_time_ratio_property(g, time_ratio):
    t_ratio = g.new_edge_property('double')
    t_ratio.a = time_ratio
    g.ep['t_ratio'] = t_ratio
    return g


def create_simple_graph(file):
    adjacency_matrix = load_npz(file).todense()
    grid_count = adjacency_matrix.shape[0]
    nnz_index = adjacency_matrix.nonzero()
    weights = adjacency_matrix[nnz_index]
    return get_base_graph(grid_count, nnz_index, weights)


def create_full_graph(adjacency_file, min_temp_file, max_temp_file):
    adjacency_matrix = load_npz(adjacency_file).todense()
    grid_count = adjacency_matrix.shape[0]
    min_temp_matrix = load_npz(min_temp_file).todense()
    max_temp_matrix = load_npz(max_temp_file).todense()

    nnz_index = adjacency_matrix.nonzero()
    print('No. of connections(nnz prob): ', len(nnz_index[0]))
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix[nnz_index]
    max_temp = max_temp_matrix[nnz_index]
    print(np.min(min_temp), np.max(max_temp))

    g = get_base_graph(grid_count, nnz_index, weights)
    g = add_temperature_properties(g, min_temp, max_temp)
    return g


def create_temp_range_graph(adjacency_file, min_temp_file, max_temp_file, temp_range):
    adjacency_matrix = load_npz(adjacency_file).todense()
    grid_count = adjacency_matrix.shape[0]
    max_temp_matrix = load_npz(max_temp_file)
    min_temp_matrix = load_npz(min_temp_file)
    temp_range_matrix = (max_temp_matrix - min_temp_matrix).todense()
    temp_filtered_matrix = np.where(temp_range_matrix > temp_range, 0, temp_range_matrix)
    filtered_matrix = np.multiply(temp_filtered_matrix, adjacency_matrix)

    nnz_index = filtered_matrix.nonzero()
    print('No. of connections(nnz prob): ', len(nnz_index[0]))
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix.todense()[nnz_index]
    max_temp = max_temp_matrix.todense()[nnz_index]

    g = get_base_graph(grid_count, nnz_index, weights)
    g = add_temperature_properties(g, min_temp, max_temp)

    # norm_min_t = g.new_edge_property('double')
    # norm_min_t.a = (min_temp - np.min(min_temp)) / np.ptp(min_temp)
    # g.ep['norm_min_t'] = norm_min_t

    return g


def create_temp_min_max_graph(adjacency_file, min_temp_file, max_temp_file, min_temp_accept, max_temp_accept):
    adjacency_matrix = load_npz(adjacency_file).todense()
    grid_count = adjacency_matrix.shape[0]
    min_temp_matrix = load_npz(min_temp_file).todense()
    max_temp_matrix = load_npz(max_temp_file).todense()

    min_temp_filter = np.where(min_temp_matrix < min_temp_accept, 0, min_temp_matrix)
    max_temp_filter = np.where(max_temp_matrix > max_temp_accept, 0, max_temp_matrix)
    t_matrix = np.multiply(min_temp_filter, max_temp_filter)
    filtered_matrix = np.multiply(t_matrix, adjacency_matrix)

    nnz_index = filtered_matrix.nonzero()
    print('No. of connections(nnz prob): ', len(nnz_index[0]))
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix[nnz_index]
    max_temp = max_temp_matrix[nnz_index]
    print(np.min(min_temp), np.max(max_temp))

    g = get_base_graph(grid_count, nnz_index, weights)
    g = add_temperature_properties(g, min_temp, max_temp)
    return g


def create_prob_filtered_graph(adjacency_file, min_temp_file, max_temp_file, prob_cutoff):
    adjacency_matrix = load_npz(adjacency_file).todense()
    grid_count = adjacency_matrix.shape[0]
    min_temp_matrix = load_npz(min_temp_file).todense()
    max_temp_matrix = load_npz(max_temp_file).todense()

    filtered_matrix = np.where(adjacency_matrix < prob_cutoff, 0, adjacency_matrix)
    # save_npz('/Users/dmanral/Desktop/Analysis/TARA/Task4/ProcessedTM/Annual_Avg_Prob_Cutoff_pt005_csr.npz',
    #          csr_matrix(filtered_matrix))

    nnz_index = filtered_matrix.nonzero()
    print('No. of connections(nnz prob): ', len(nnz_index[0]))
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix[nnz_index]
    max_temp = max_temp_matrix[nnz_index]
    print(np.min(min_temp), np.max(max_temp))

    g = get_base_graph(grid_count, nnz_index, weights)
    g = add_temperature_properties(g, min_temp, max_temp)
    return g


def create_minimum_time_graph(adjacency_file, t_ratio_file):
    adjacency_matrix = load_npz(adjacency_file).todense()
    grid_count = adjacency_matrix.shape[0]
    nnz_index = adjacency_matrix.nonzero()
    weights = adjacency_matrix[nnz_index]
    g = get_base_graph(grid_count, nnz_index, weights)
    time_ratio = load_npz(t_ratio_file).todense()

    g = add_minimum_time_ratio_property(g, time_ratio[nnz_index])
    return g


def get_most_probable_path(g, s, d):
    vlist, elist = gt.shortest_path(g, s, d, weights=g.ep['weight'])
    path = [int(v) for v in vlist]
    # print(path)
    # if path:
    #     get_path_probabilities(g, path)
    #     try:
    #         temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
    #         print(np.around(temp_range, 2))
    #     except KeyError:
    #         print("Temperature not included in the graph")
    #     print('path length:', len(path))
    return path


def get_path_probabilities(g, path):
    probs = np.zeros((len(path) - 1))
    for i in range(len(path) - 1):
        v0 = path[i]
        v1 = path[i + 1]
        prob0 = get_prob(g, v0, v1)
        probs[i] = np.round(prob0, 4)
    print(probs)


def get_shortest_path(g, s, d):
    vlist, elist = gt.shortest_path(g, s, d)
    path = [int(v) for v in vlist]
    # if path:
    #     print(len(path))
    #     print(path)
    #     get_path_probabilities(g, path)
    #     try:
    #         temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
    #         print(np.around(temp_range, 2))
    #     except KeyError:
    #         print("Temperature not included in the graph")
    return path


def get_shortest_paths_subset(g, s, d, path_length):
    cnt = gt.count_shortest_paths(g, s, d)
    print('Total paths: ', cnt)
    if cnt > 100:
        count = 100
    else:
        count = cnt
    paths = np.empty((count, path_length), dtype=np.int32)
    for i in range(count):
        paths[i] = gt.random_shortest_path(g, s, d)
    print('paths computed. Count: ', count)
    return paths


def get_all_shortest_paths(g, s, d, path_length):
    count = gt.count_shortest_paths(g, s, d)
    print('Total paths: ', count)
    paths = gt.all_shortest_paths(g, s, d)
    all_paths = np.empty((count, path_length), dtype=np.int32)
    for i, p in zip(range(count), paths):
        all_paths[i] = [int(v) for v in p]
    return all_paths


def get_prob(g, s, d):
    e = g.edge(s, d)
    if e is not None:
        e_idx = g.edge_index[e]
        return g.ep['probability'].a[e_idx]
    return 0


def get_time_from_most_probable_path(g, path):
    t = 0
    time_laps = np.empty(len(path))
    time_laps[0] = 0

    for i in range(len(path) - 1):
        v0 = path[i]
        v1 = path[i + 1]
        prob0 = get_prob(g, v0, v0)
        prob1 = get_prob(g, v0, v1)
        t += (prob0 / prob1) + 1
        # print(np.round(prob0, 4), np.round(prob1, 4), np.round(t / 12, 4))
        time_laps[i + 1] = t / 12

    # print('Total time in years:', t / 12)
    return time_laps


# def get_least_probable_path(g, s, d):
#     vlist, elist = gt.shortest_path(g, s, d, weights=g.ep['probability'])
#     path = [int(v) for v in vlist]
#     print(path)
#     get_path_probabilities(g, path)
#     try:
#         temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
#         print(np.around(temp_range, 2))
#     except KeyError:
#         print("Temperature not included in the graph")
#     print(len(path))
#     return path


def get_all_paths(g, s, d, path_length):
    count = gt.count_shortest_paths(g, s, d)
    all_paths = gt.all_shortest_paths(g, s, d)
    # gt.shortest_distance(g, s, d, weights=g.ep['weight'], directed=True)
    paths = np.empty((count, path_length), dtype=np.int32)
    for i, p in range(count), all_paths:
        paths[i] = p
    print('paths computed. Count: ', count)
    return paths


# def get_time_from_shortest_path(g, path):
#     t = 0
#     time_laps = np.empty(len(path) - 1)
#
#     for i in range(len(path) - 1):
#         v0 = path[i]
#         v1 = path[i + 1]
#         prob0 = get_prob(g, v0, v0)
#         prob1 = get_prob(g, v0, v1)
#         t += prob0 / (prob1 + 1)
#         print(np.round(prob0, 4), np.round(prob1, 4), np.round(t / 12, 4))
#         time_laps[i] = t / 12
#
#     print('Total time in years:', t / 12)
#     return time_laps
#

# def get_min_time_of_all_paths(g, paths):
#
def check_if_edge_exists(g, s, d):
    ed = g.edge(s, d)
    if ed:
        return True
    return False


def get_shortest_path_min_temp_range(g, s, d):
    vlist, elist = gt.shortest_path(g, s, d, weights=g.ep['norm_min_t'])  # Bellman-Ford algorithm.
    path = [int(v) for v in vlist]
    if path:
        print(len(path))
        print(path)
        get_path_probabilities(g, path)
        try:
            temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
            print(np.around(temp_range, 2))
        except KeyError:
            print("Temperature not included in the graph")
    return path


def incoming_connections(g, idx):
    eds = g.get_in_edges(idx, [g.ep['probability'], g.ep['min_t'], g.ep['max_t']])
    return eds[:, 0].astype(np.int32), eds[:, 2], eds[:, 3], eds[:, 4]


def outgoing_connections(g, idx):
    eds = g.get_out_edges(idx, [g.ep['probability'], g.ep['min_t'], g.ep['max_t']])
    return eds[:, 1].astype(np.int32), eds[:, 2], eds[:, 3], eds[:, 4]


def minimum_time_path_malley2021(g, s, d):
    vlist, elist = gt.shortest_path(g, s, d, weights=g.ep['t_ratio'])
    path = [int(v) for v in vlist]
    # if path:
    #     print(len(path))
    #     print(path)
    #     get_path_probabilities(g, path)
    #     try:
    #         temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
    #         print(np.around(temp_range, 2))
    #     except KeyError:
    #         print("Temperature not included in the graph")
    return path
