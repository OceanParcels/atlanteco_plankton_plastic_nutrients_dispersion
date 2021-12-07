import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task6_Sens/Resolution4/'


def get_empty_array(mons, no_grids):
    return np.zeros((mons, no_grids, no_grids + 2))


def get_dense_matrix(matrix_data, indices, indptr, no_grids):
    return csr_matrix((matrix_data, indices, indptr), shape=(no_grids, no_grids + 2)).todense()


def compute_grid_annual_average(matrix, no_grids):
    avg_matrix = np.mean(matrix, axis=0)
    print(np.min(avg_matrix), np.max(avg_matrix))
    print(np.min(avg_matrix[:, :no_grids]), np.max(avg_matrix[:, :no_grids]))
    assert np.sum(avg_matrix, axis=1).all() == 1
    save_npz(home_folder + 'Annual_Avg_FullAdjacency_csr.npz', csr_matrix(avg_matrix))
    save_npz(home_folder + 'Annual_Avg_DomainAdjacency_csr.npz', csr_matrix(avg_matrix[:, :no_grids]))


def compute_grid_nnz_average(matrix, nnz_count_matrix, matrix_type):
    sum_matrix = np.sum(matrix, axis=0)
    avg_matrix = sum_matrix / nnz_count_matrix
    print(matrix_type, np.nanmin(avg_matrix), np.nanmax(avg_matrix))
    save_npz(home_folder + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))


def compute_total_transitions_field_average(matrix, transition_count, matrix_type):
    sum_matrix = np.sum(matrix, axis=0)
    avg_matrix = sum_matrix / transition_count
    print(matrix_type, np.nanmin(avg_matrix), np.nanmax(avg_matrix))
    save_npz(home_folder + 'Annual_SUM_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))


def main():
    # no_grids = 8245
    no_grids = len(pd.read_csv(home_folder + 'MasterHexList_Res4.csv',header=None))
    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    month_count = len(months)

    monthly_full_trans_prob = get_empty_array(month_count, no_grids)
    monthly_full_min_temp = get_empty_array(month_count, no_grids)
    monthly_full_max_temp = get_empty_array(month_count, no_grids)
    monthly_full_min_sal = get_empty_array(month_count, no_grids)
    monthly_full_max_sal = get_empty_array(month_count, no_grids)

    # Saving format: from MonthlySparseTM.py
    # np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
    #                         mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
    #                         maxsal=avg_max_sal_per_grid, indices=mon_trans_matrix.indices,
    #                         indptr=mon_trans_matrix.indptr)

    for mon, i in zip(months, range(len(months))):
        with np.load(home_folder + 'CSR_{0}.npz'.format(mon)) as data:
            monthly_full_trans_prob[i] = get_dense_matrix(data['transprob'], data['indices'], data['indptr'], no_grids)
            monthly_full_min_temp[i] = get_dense_matrix(data['mintemp'], data['indices'], data['indptr'], no_grids)
            monthly_full_max_temp[i] = get_dense_matrix(data['maxtemp'], data['indices'], data['indptr'], no_grids)
            monthly_full_min_sal[i] = get_dense_matrix(data['minsal'], data['indices'], data['indptr'], no_grids)
            monthly_full_max_sal[i] = get_dense_matrix(data['maxsal'], data['indices'], data['indptr'], no_grids)

    # Compute normalized avg or sum of transitions
    option = 1
    if option == 1:
        compute_grid_annual_average(monthly_full_trans_prob, no_grids)
        # remaining without the new column(8245) and deleted(8246)
        nnz_count_matrix = np.count_nonzero(monthly_full_trans_prob[:, :, :no_grids], axis=0)
        nnz_count_matrix[nnz_count_matrix == 0] = 1
        compute_grid_nnz_average(monthly_full_min_temp[:, :, :no_grids], nnz_count_matrix, 'MinTemperature')
        compute_grid_nnz_average(monthly_full_max_temp[:, :, :no_grids], nnz_count_matrix, 'MaxTemperature')
        compute_grid_nnz_average(monthly_full_min_sal[:, :, :no_grids], nnz_count_matrix, 'MinSalinity')
        compute_grid_nnz_average(monthly_full_max_sal[:, :, :no_grids], nnz_count_matrix, 'MaxSalinity')
    elif option == 2:
        annual_sum_matrix = np.sum(monthly_full_trans_prob, axis=0)
        print("Min/Max Total transitions", np.min(annual_sum_matrix), np.max(annual_sum_matrix))
        save_npz(home_folder + 'Annual_SUM_FullAdjacency_csr.npz', csr_matrix(annual_sum_matrix))
        save_npz(home_folder + 'Annual_SUM_DomainAdjacency_csr.npz', csr_matrix(annual_sum_matrix[:, :no_grids]))

        nnz_count_matrix = annual_sum_matrix[:, :no_grids]
        nnz_count_matrix[nnz_count_matrix == 0] = 1

        compute_total_transitions_field_average(monthly_full_min_temp[:, :, :no_grids], nnz_count_matrix,
                                                'MinTemperature')
        compute_total_transitions_field_average(monthly_full_max_temp[:, :, :no_grids], nnz_count_matrix,
                                                'MaxTemperature')
        compute_total_transitions_field_average(monthly_full_min_sal[:, :, :no_grids], nnz_count_matrix,
                                                'MinSalinity')
        compute_total_transitions_field_average(monthly_full_max_sal[:, :, :no_grids], nnz_count_matrix,
                                                'MaxSalinity')
    else:
        raise ValueError('option value is incorrect')


if __name__ == '__main__':
    main()