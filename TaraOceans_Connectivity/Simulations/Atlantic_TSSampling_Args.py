import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode
import numpy as np
from kernels.samplefields import SampleTSFields
import sys

    
def delete_particle(particle, fieldset, time):
    particle.delete()
    
# verify input parameters: arg1= Year, arg2= Month
args=sys.argv
assert len(args) == 3   

start_year = np.int32(args[1])
assert 2009 <= start_year <= 2018

start_mon = np.int32(args[2]) 
assert 1 <= start_mon <= 12

if start_mon == 12:
    if start_year == 2018:
        raise ValueError('data unavailable for complete simulation.')
    end_mon = 1
    end_year = start_year + 1
else:
    end_mon = start_mon + 1
    end_year = start_year


data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, 1, 12, 0, 0)
simulation_end = datetime(end_year, end_mon, 1, 12, 0, 0)
print(simulation_start, simulation_end)

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

mon = simulation_start.strftime("%b")
output_file = pset.ParticleFile(name="/scratch/manra003/tara_res5_01/FullTara_Res5_TS_{0}{1}_dt600.nc".format(mon, start_year), outputdt=timedelta((simulation_end-simulation_start).days))

sample_kernel = pset.Kernel(SampleTSFields)

pset.execute(AdvectionRK4 + sample_kernel,                
             endtime=simulation_end,
             dt=600,                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()