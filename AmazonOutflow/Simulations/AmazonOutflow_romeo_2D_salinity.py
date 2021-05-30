import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import numpy as np

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

year = 2009

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}0[7-9]*_U.nc'.format(year)) + glob(data_path + 'ROMEO.01_1d_uo_{0}10*_U.nc'.format(year)))
vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}0[7-9]*_V.nc'.format(year)) + glob(data_path + 'ROMEO.01_1d_vo_{0}10*_V.nc'.format(year)))
sfiles =  sorted(glob(data_path + 'ROMEO.01_1d_so_{0}0[7-9]*_T.nc'.format(year)) + glob(data_path + 'ROMEO.01_1d_so_{0}10*_T.nc'.format(year)))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles},
             'S': {'lon': mesh_mask, 'lat': mesh_mask, 'data': sfiles}}

variables = {'U': 'uo',
             'V': 'vo',
             'S': 'so'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

simulation_start = datetime(year, 7, 1, 12, 0, 0)
simulation_end = datetime(year, 10, 31, 12, 0, 0)

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')


class Particle(JITParticle):
    s = Variable('s', dtype=np.float32)

def SampleSalinity(particle, fieldset, time):
    particle.s = fieldset.S[time, particle.depth, particle.lat, particle.lon]
    
pset = ParticleSet.from_line(fieldset=fieldset, 
                             size=200, 
                             pclass=Particle,
                             start=(-49.0664412, 0.930062), 
                             finish=(-48.8676954, 0.6292335), 
                             repeatdt=timedelta(days=1),
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/amazon/Atlanteco_Romeo_July-Oct{0}_OFF50km_z0m_sal.nc".format(year), outputdt=timedelta(days=1))

Salinity_kernel = pset.Kernel(SampleSalinity)

pset.execute(AdvectionRK4 + Salinity_kernel,                
             endtime=simulation_end,
             dt=300,                       
             output_file=output_file)

output_file.close()

