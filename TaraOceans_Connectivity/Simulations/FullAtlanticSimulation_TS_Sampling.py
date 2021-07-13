import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode
import numpy as np

  
def SampleFields(particle, fieldset, time):
    temp = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    sal = fieldset.S[time, particle.depth, particle.lat, particle.lon]

    if time == 0.0:
        particle.min_temp = temp
        particle.max_temp = temp
        particle.min_sal = sal
        particle.max_sal = sal
        
    if temp < particle.min_temp:
        particle.min_temp = temp
    elif temp > particle.max_temp:
        particle.max_temp = temp
    
    if sal < particle.min_sal:
        particle.min_sal = sal
    elif sal > particle.max_sal:
        particle.max_sal = sal
    
def delete_particle(particle, fieldset, time):
    particle.delete()
    

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(2011, 1, 1, 12, 0, 0)
simulation_end = datetime(2011, 2, 1, 12, 0, 0)

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_uo_{0}{1}01_grid_U.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_vo_{0}{1}01_grid_V.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

tfiles =  sorted(glob(data_path + 'ROMEO.01_1d_thetao_{0}{1}*_T.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_thetao_{0}{1}01_grid_T.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

sfiles =  sorted(glob(data_path + 'ROMEO.01_1d_so_{0}{1}*_T.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))) + \
                 glob(data_path + 'ROMEO.01_1d_so_{0}{1}01_grid_T.nc'.\
                      format(simulation_end.strftime("%Y"), simulation_end.strftime("%m"))))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles},
             'T': {'lon': mesh_mask, 'lat': mesh_mask, 'data': tfiles},
             'S': {'lon': mesh_mask, 'lat': mesh_mask, 'data': sfiles}
            }

variables = {'U': 'uo',
             'V': 'vo',
             'T': 'thetao',
             'S': 'so'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize=False)

coords = pd.read_csv(r'/scratch/manra003/data/Nemo_H3Release_LatLon_Res5.csv')

class Particle(JITParticle):
    min_temp = Variable('min_temp', dtype=np.float32)
    max_temp = Variable('max_temp', dtype=np.float32) 
    min_sal = Variable('min_sal', dtype=np.float32)
    max_sal = Variable('max_sal', dtype=np.float32)  
    
    
pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=Particle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/tara_res5_01/FullTara_Res5_TS_Jan2011_dt600.nc", outputdt=timedelta((simulation_end-simulation_start).days))

sample_kernel = pset.Kernel(SampleFields)

pset.execute(AdvectionRK4 + sample_kernel,                
             endtime=simulation_end,
             dt=600,                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()