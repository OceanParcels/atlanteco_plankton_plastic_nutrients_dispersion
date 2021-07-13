import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime
import pandas as pd

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(2009, 9, 1, 12, 0, 0)
simulation_end = datetime(2009, 9, 16, 12, 0, 0)

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

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

coords = pd.read_csv(r'/scratch/manra003/data/Amazon_Retro_0625_FULL.csv')

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/Amazon_Retro_0625_FULL_Sep01-16_2009_odt6h.nc", outputdt=timedelta(hours=6))

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=300,                       
             output_file=output_file)

output_file.close()

