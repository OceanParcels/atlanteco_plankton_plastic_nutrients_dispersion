from graph_tool.all import Graph, shortest_path, count_shortest_paths, random_shortest_path
from scipy.sparse import load_npz
import numpy as np


def get_base_graph(nnz_index, weights):
    g = Graph()
    g.add_edge_list(np.transpose(nnz_index))
    ew = g.new_edge_property('double')
    ew.a = -np.log(weights)
    g.ep['weight'] = ew

    eprob = g.new_edge_property('double')
    eprob.a = weights
    g.ep['probability'] = eprob
    return g


def create_simple_graph(file):
    adjacency_matrix = load_npz(file).todense()
    nnz_index = adjacency_matrix.nonzero()
    weights = adjacency_matrix[nnz_index]
    return get_base_graph(nnz_index, weights)


def create_temp_range_graph(adjacency_file, min_temp_file, max_temp_file, temp_range):
    adjacency_matrix = load_npz(adjacency_file).todense()
    max_temp_matrix = load_npz(max_temp_file)
    min_temp_matrix = load_npz(min_temp_file)
    temp_range_matrix = (max_temp_matrix - min_temp_matrix).todense()
    filtered_matrix = np.where(temp_range_matrix > temp_range, 0, temp_range_matrix)

    nnz_index = filtered_matrix.nonzero()
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix.todense()[nnz_index]
    max_temp = max_temp_matrix.todense()[nnz_index]

    g = get_base_graph(nnz_index, weights)
    min_t = g.new_edge_property('double')
    min_t.a = min_temp
    g.ep['min_t'] = min_t
    max_t = g.new_edge_property('double')
    max_t.a = max_temp
    g.ep['max_t'] = max_t
    return g


def create_temp_min_max_graph(adjacency_file, min_temp_file, max_temp_file, min_temp_accept, max_temp_accept):
    adjacency_matrix = load_npz(adjacency_file).todense()
    min_temp_matrix = load_npz(min_temp_file).todense()
    max_temp_matrix = load_npz(max_temp_file).todense()

    min_temp_filter = np.where(min_temp_matrix < min_temp_accept, 0, min_temp_matrix)
    max_temp_filter = np.where(max_temp_matrix > max_temp_accept, 0, max_temp_matrix)
    filtered_matrix = np.multiply(min_temp_filter, max_temp_filter)

    nnz_index = filtered_matrix.nonzero()
    weights = adjacency_matrix[nnz_index]
    min_temp = min_temp_matrix[nnz_index]
    max_temp = max_temp_matrix[nnz_index]
    print(np.min(min_temp), np.max(max_temp))

    g = get_base_graph(nnz_index, weights)
    min_t = g.new_edge_property('double')
    min_t.a = min_temp
    g.ep['min_t'] = min_t
    max_t = g.new_edge_property('double')
    max_t.a = max_temp
    g.ep['max_t'] = max_t
    return g


def get_most_probable_path(g, s, d):
    vlist, elist = shortest_path(g, s, d, weights=g.ep['weight'])
    path = [int(v) for v in vlist]
    print(path)
    try:
        temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
        print(np.around(temp_range, 2))
    except KeyError:
        print("Temperature not included in the graph")
    print(len(path))
    return path


def get_shortest_path(g, s, d):
    vlist, elist = shortest_path(g, s, d)
    path = [int(v) for v in vlist]
    print(path)
    try:
        temp_range = [(g.ep['min_t'][e], g.ep['max_t'][e]) for e in elist]
        print(np.around(temp_range, 2))
    except KeyError:
        print("Temperature not included in the graph")
    print(len(path))
    return path


def get_shortest_paths_subset(g, s, d, path_length):
    cnt = count_shortest_paths(g, s, d)
    print('Total paths: ', cnt)
    if cnt > 100:
        count = 100
    else:
        count = cnt
    paths = np.empty((count, path_length), dtype=np.int32)
    for i in range(count):
        paths[i] = random_shortest_path(g, s, d)
    print('paths computed. Count: ', count)
    return paths


def get_time_from_most_probable_path(g, path):
    t = 0
    time_laps = np.empty(len(path) - 1)

    def get_prob(s, d):
        e = g.edge(s, d)
        if e is not None:
            e_ind = g.edge_index[e]
            return g.ep['probability'].a[e_ind]
        return 0

    for i in range(len(path) - 1):
        v0 = path[i]
        v1 = path[i + 1]
        prob0 = get_prob(v0, v0)
        prob1 = get_prob(v0, v1)
        t += (prob0 / prob1) + 1
        print(np.round(prob0, 4), np.round(prob1, 4), np.round(t / 12, 4))
        time_laps[i] = t / 12

    print('Total time in years:', t / 12)
    return time_laps
