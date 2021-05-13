import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime
import pandas as pd
import random
from parcels.tools.statuscodes import ErrorCode

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(2010, 12, 15, 12, 0, 0)
simulation_end = datetime(2011, 1, 15, 12, 0, 0)

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_uo_{0}{1}15_grid_U.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_vo_{0}{1}15_grid_V.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')

seed_points = pd.read_csv(r'/scratch/manra003/Nemo_H3Release_LatLon.csv')

# Randomly select a lat-lon from each hexagon/pentagon: Total: 8114
hex_list = seed_points["Res3_HexId"]
unique_hex_list = hex_list.unique()

seed_pool = seed_points.groupby(["Res3_HexId"]).indices

indexes = [random.choice(seed_pool[hex]) for hex in unique_hex_list]

seed_lats = [seed_points['Latitudes'].loc[i] for i in indexes]
seed_lons = [seed_points['Longitudes'].loc[i] for i in indexes]

# # region: standard release points for sensitivity analysis
# seed_points = pd.read_csv(r'/scratch/manra003/Analysis_Sample_Nemo_H3Release_LatLon.csv')
# seed_lats = seed_points['Latitudes']
# seed_lons = seed_points['Longitudes']
# # endregion

assert len(seed_lats) == 8114

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=seed_lons,
                             lat=seed_lats,
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/tara_data_15/FullAtlantic_2D_15Dec2010_1month.nc", outputdt=timedelta((simulation_end-simulation_start).days))

def delete_particle(particle, fieldset, time):
    # track particles that are removed
    print("id: %d, lat: %f, lon: %f, time: %f" % (particle.id, particle.lat, particle.lon, time))    
    particle.delete()

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=300,                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()

