import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Variable
from datetime import timedelta, datetime
import numpy as np

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_20180[7-9]*_U.nc') + glob(data_path + 'ROMEO.01_1d_uo_20181[0-1]*_U.nc'))
vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_20180[7-9]*_V.nc') + glob(data_path + 'ROMEO.01_1d_vo_20181[0-1]*_V.nc'))
sfiles =  sorted(glob(data_path + 'ROMEO.01_1d_so_20180[7-9]*_T.nc') + glob(data_path + 'ROMEO.01_1d_so_20181[0-1]*_T.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles},
             'S': {'lon': mesh_mask, 'lat': mesh_mask, 'data': sfiles}}

variables = {'U': 'uo',
             'V': 'vo',
             'S': 'so'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

simulation_start = datetime(2018, 7, 1, 12, 0, 0)
simulation_end = datetime(2018, 11, 30, 12, 0, 0)

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
                             start=(-46.3224637, 2.5242153), 
                             finish=(-46.1222564, 2.2223863), 
                             repeatdt=timedelta(days=1),
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/Atlanteco_Romeo_July-Nov2018_OFF400km_z0m_sal.nc", outputdt=timedelta(days=1))

Salinity_kernel = pset.Kernel(SampleSalinity)

pset.execute(AdvectionRK4 + Salinity_kernel,                
             endtime=simulation_end,
             dt=300,                       
             output_file=output_file)

output_file.close()

