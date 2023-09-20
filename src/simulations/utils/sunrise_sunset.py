import numpy as np
import netCDF4 as nc
import xarray as xr
import pandas as pd
from datetime import date
import calendar

def load_sunrise_sunset(days, folder_path):
    '''
    Function to load the sunrise and sunset files computed using SunriseSunset Calculator tool.
    The data is then replicated for the years for which simulation needs to run.
    Hence requires the 'dates' used to load the ocean model data
    '''
    years= np.unique([d.year for d in days])

    # load a non leap-year values for sunrise and sunset
    sunrise_nc = nc.Dataset(folder_path + 'SunriseTime_2x2_1d_2015.nc','r')
    sunrise_da = xr.DataArray(sunrise_nc['sunrise'][::], coords={"time": sunrise_nc['time'][:], "lat": sunrise_nc['lat'][:], "lon": sunrise_nc['lon'][:]})

    sunset_nc = nc.Dataset(folder_path + 'SunsetTime_2x2_1d_2015.nc','r')
    sunset_da = xr.DataArray(sunset_nc['sunset'][::], coords={"time": sunset_nc['time'][:], "lat": sunset_nc['lat'][:], "lon": sunset_nc['lon'][:]})
    
    # for unique years in the data files, get the datetimes for all days and remove the leap year dates. 
    # It is okay to have missing values for a day in Feb, interpolation should work fine. It is just o have consistent values for all days
    #remove leap year dates from the list
    dates = pd.date_range(date(years[0],1,1), date(years[0],12,31))
    if calendar.isleap(years[0]):
        dates = dates.delete(31+28)
    sunrise_da['time']= dates
    sunset_da['time']= dates
    sunrise_nc.close()
    sunset_nc.close()
    sunrise_copy=sunrise_da.copy()
    sunset_copy=sunset_da.copy()

    for y in years[1:]:
        dates = pd.date_range(date(y, 1, 1), date(y, 12, 31))
        if calendar.isleap(y):
            dates = dates.delete(31+28)
        sunrise_copy['time']= dates
        sunrise_da = xr.combine_by_coords([sunrise_da, sunrise_copy], coords=['time'])

        sunset_copy['time']= dates
        sunset_da = xr.combine_by_coords([sunset_da, sunset_copy], coords=['time'])

    return sunrise_da, sunset_da
