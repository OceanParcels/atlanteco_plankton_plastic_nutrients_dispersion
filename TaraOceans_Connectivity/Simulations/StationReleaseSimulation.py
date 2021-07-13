import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
from parcels.tools.statuscodes import ErrorCode
import numpy as np
import pandas as pd

 
# def SampleTemperature(particle, fieldset, time):
#     particle.t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    
    
def delete_particle(particle, fieldset, time):
#     print("id: %d, lat: %f, lon: %f, time: %f" % (particle.id, particle.lat, particle.lon, time)) 
    particle.delete()


data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_201[0-5]*_U.nc'))
vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_201[0-5]*_V.nc'))
# tfiles =  sorted(glob(data_path + 'ROMEO.01_1d_thetao_*_T.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}
#              'T': {'lon': mesh_mask, 'lat': mesh_mask, 'data': tfiles}}

variables = {'U': 'uo',
             'V': 'vo'}
#              'T': 'thetao'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

simulation_start = datetime(2010, 1, 1, 12, 0, 0)
simulation_end = datetime(2015, 12 , 31, 12, 0, 0)

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')


# class Particle(JITParticle):
#     t = Variable('t', dtype=np.float32)
    

coords = pd.read_csv(r'/scratch/manra003/data/H3_Res9points_OA236.csv')

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start)
 
output_file = pset.ParticleFile(name="/scratch/manra003/Tara_OA236_Res9_6years_dt10mins_out3h.nc", outputdt=timedelta(hours=3))
print('kernel started')
# Temperature_kernel = pset.Kernel(SampleTemperature)
    
pset.execute(AdvectionRK4, # + Temperature_kernel,                
             endtime=simulation_end,
             dt=600,                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()