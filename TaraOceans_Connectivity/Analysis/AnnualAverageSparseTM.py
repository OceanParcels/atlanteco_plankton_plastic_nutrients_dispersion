import numpy as np
from scipy.sparse import csr_matrix, save_npz

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/ProcessedTM/'
no_particles = 377583
no_grids = 8245


def main():
    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    def get_empty_array():
        return np.zeros((len(months), no_grids, no_grids + 2))

    monthly_full_trans_prob = get_empty_array()
    monthly_full_min_temp = get_empty_array()
    monthly_full_max_temp = get_empty_array()
    monthly_full_min_sal = get_empty_array()
    monthly_full_max_sal = get_empty_array()

    # Saving format: from MonthlySparseTM.py
    # np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
    #                         mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
    #                         maxsal=avg_max_sal_per_grid, indices=mon_trans_matrix.indices,
    #                         indptr=mon_trans_matrix.indptr)

    def get_dense_matrix(matrix_data):
        return csr_matrix((matrix_data, data['indices'], data['indptr']),
                          shape=(no_grids, no_grids + 2)).todense()

    for mon, i in zip(months, range(len(months))):
        with np.load(home_folder + 'CSR_{0}.npz'.format(mon)) as data:
            monthly_full_trans_prob[i] = get_dense_matrix(data['transprob'])
            monthly_full_min_temp[i] = get_dense_matrix(data['mintemp'])
            monthly_full_max_temp[i] = get_dense_matrix(data['maxtemp'])
            monthly_full_min_sal[i] = get_dense_matrix(data['minsal'])
            monthly_full_max_sal[i] = get_dense_matrix(data['maxsal'])

    def compute_grid_annual_average(matrix, matrix_type):
        avg_matrix = np.average(matrix, axis=0)
        print(matrix_type, np.min(avg_matrix), np.max(avg_matrix))
        save_npz(home_folder + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))

    compute_grid_annual_average(monthly_full_trans_prob, 'FullAdjacency')
    # remaining without the new column(8245) and deleted(8246)
    compute_grid_annual_average(monthly_full_trans_prob[:, :, :no_grids], 'DomainAdjacency')
    compute_grid_annual_average(monthly_full_min_temp[:, :, :no_grids], 'MinTemperature')
    compute_grid_annual_average(monthly_full_max_temp[:, :, :no_grids], 'MaxTemperature')
    compute_grid_annual_average(monthly_full_min_sal[:, :, :no_grids], 'MinSalinity')
    compute_grid_annual_average(monthly_full_max_sal[:, :, :no_grids], 'MaxSalinity')


if __name__ == '__main__':
    main()
