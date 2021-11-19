import pandas as pd
import re
import h3


def convert_format(deg_min):
    """
    convert deg_min to deg_decimal lat/lon
    """
    dir = -1 if re.search('[swSW]', deg_min) else 1
    out = re.split("[°′\"]", deg_min)
    minute = out[1] if len(out) > 2 else 0
    second = out[2] if len(out) > 3 else 0

    return dir * (float(out[0]) + float(minute) / 60 + float(second) / 3600)


home_folder = '/Users/dmanral/Desktop/Analysis/UvA/'
stations_data = pd.read_excel(home_folder + 'Sampling.xls', header=0)

org_lats = stations_data['Latitude']
org_lons = stations_data['Longitude']

lats = [convert_format(i) for i in org_lats]
lons = [convert_format(i) for i in org_lons]

# get hex mapping for stations
hex_list = [h3.geo_to_h3(x, y, 3) for x, y in zip(lats, lons)]

df = pd.DataFrame({'Code': stations_data['Station'],
                   'Latitude': lats,
                   'Longitude': lons,
                   'H3Id': hex_list})
df.to_csv(home_folder + 'StationsUvA.csv')
