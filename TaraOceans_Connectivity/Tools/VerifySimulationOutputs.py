import numpy as np
from glob import glob
import xarray as xr
import pandas as pd

home_folder = '/scratch/dmanral/'
seed_points = pd.read_csv(home_folder + 'data/Nemo_H3Release_LatLon_Res5.csv')
master_hexId = seed_points['Res3_HexId'].unique()
master_particleId = np.arange(0, len(master_hexId), 1)
print(len(seed_points['Latitudes']),len(master_hexId))

run_for_months = 12
total_no_particles = len(seed_points['Latitudes'])

months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
years = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]

def verify_date(date, day, month, yr):
    d = pd.to_datetime(date)
    if d.day == day and d.month == month and d.year == yr:
        return date
    raise ValueError
    
for mon in months[:]:
    # files from all the years for that month
    monthly_files = sorted(glob(home_folder + 'tara_res5_01/FullTara_Res5_TS_{0}*_dt600.nc'.format(mon)))
    mon_index = np.where(months == mon)[0]
    print(len(monthly_files))
    for file, year in zip(monthly_files, years):
        ds = xr.open_dataset(file).load()

        # ensure start date and end date are as expected
        t_min = verify_date(ds['time'][5000,0].values, 1, mon_index + 1, year)
        if mon_index < 11:
            t_max = verify_date(ds['time'][5000,1].values, 1, mon_index + 2, year)
        else:
            if year != 2018:
                t_max = verify_date(ds['time'][5000,1].values, 1, 1, year + 1)
            else:
                t_max = verify_date(ds['time'][5000,1].values, 31, mon_index + 1, year)
        print(ds['time'][5000,0].values, ds['time'][5000,1].values)
        ds.close()

print("completed")    