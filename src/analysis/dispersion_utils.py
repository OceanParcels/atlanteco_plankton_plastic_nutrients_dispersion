import numpy as np
from scipy import stats


def distance(lon1, lat1, lon2, lat2, r=6378):
    """
    Parallelised code from - Michael Denes
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    # if np.isnan([lon1, lon2, lat1, lat2]).any():
    #     return mask_value
    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    dlat = np.subtract(lat2, lat1)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2),  
                          np.multiply(np.cos(lat1), 
                                      np.multiply(np.cos(lat2), 
                                                  np.power(np.sin(np.divide(dlon, 2)), 2))))
    c = np.multiply(2, np.arcsin(np.sqrt(a)))

    return c*r


def threshold_days(array, delta_d, mask_value=104.0):
    def delta_value_index(row, tc):
        search = np.where(row>=tc)[0]
        if search.size>0:
            return search[0]
        return mask_value # masking value-  not crossed in 100 days value just one bin higher 

    days=np.empty((array.shape[0]))
    for i in range(array.shape[0]):
        days[i]=delta_value_index(array[i,:], delta_d)
    return days


def get_diff_CDF_PDF(array, delta_d, bins):
    diff = threshold_days(array, delta_d)
    print(stats.describe(diff, axis= None)) #, nan_policy='omit'
    count, _ = np.histogram(diff, bins=bins)

    pdf = count/np.sum(count)
    cdf = np.cumsum(pdf)
    return cdf, pdf