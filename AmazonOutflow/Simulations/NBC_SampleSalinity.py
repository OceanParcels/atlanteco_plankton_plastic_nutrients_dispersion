import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime
import pandas as pd

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(2009, 8, 1, 12, 0, 0)

sfiles =  data_path + 'ROMEO.01_1d_so_{0}{1}*_T.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))

filenames = {'S': {'lon': mesh_mask, 'lat': mesh_mask, 'data': sfiles}}

variables = {'S': 'so'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')

coords = pd.read_csv(r'/scratch/manra003/data/Amazon_Retro_0625_FULL.csv')

class Particle(JITParticle):
    s = Variable('s', dtype=np.float32)

def SampleSalinity(particle, fieldset, time):
    particle.s = fieldset.S[time, particle.depth, particle.lat, particle.lon]

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/Amazon_Retro_0625_salinty_Aug01_2009.nc", outputdt=timedelta(hours=1))

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=timedelta(hours=1),                       
             output_file=output_file)

output_file.close()


