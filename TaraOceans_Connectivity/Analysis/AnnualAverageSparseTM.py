import numpy as np
from scipy.sparse import csr_matrix, save_npz

home_folder = '/TARA/Task4/ProcessedTM/'


def get_empty_array(mons, no_grids):
    return np.zeros((mons, no_grids, no_grids + 2))


def get_dense_matrix(matrix_data, indices, indptr, no_grids):
    return csr_matrix((matrix_data, indices, indptr), shape=(no_grids, no_grids + 2)).todense()


def compute_grid_annual_average(matrix, matrix_type):
    avg_matrix = np.average(matrix, axis=0)
    print(matrix_type, np.min(avg_matrix), np.max(avg_matrix))
    save_npz(home_folder + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))
    

def compute_grid_nnz_average(matrix, nnz_count_matrix, matrix_type):
    sum_matrix = np.sum(matrix, axis=0)
    avg_matrix = sum_matrix / nnz_count_matrix
    print(matrix_type, np.nanmin(avg_matrix), np.nanmax(avg_matrix))
    save_npz(home_folder + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))


def main():
    no_grids = 8245
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

    compute_grid_annual_average(monthly_full_trans_prob, 'FullAdjacency')
    # remaining without the new column(8245) and deleted(8246)
    compute_grid_annual_average(monthly_full_trans_prob[:, :, :no_grids], 'DomainAdjacency')
    nnz_count_matrix = np.count_nonzero(monthly_full_trans_prob[:, :, :no_grids], axis=0)
    nnz_count_matrix[nnz_count_matrix == 0] = 1
    compute_grid_nnz_average(monthly_full_min_temp[:, :, :no_grids], nnz_count_matrix, 'MinTemperature')
    compute_grid_nnz_average(monthly_full_max_temp[:, :, :no_grids], nnz_count_matrix, 'MaxTemperature')
    compute_grid_nnz_average(monthly_full_min_sal[:, :, :no_grids], nnz_count_matrix, 'MinSalinity')
    compute_grid_nnz_average(monthly_full_max_sal[:, :, :no_grids], nnz_count_matrix, 'MaxSalinity')


if __name__ == '__main__':
    main()