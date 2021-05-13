import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_201807*_U.nc') + glob(data_path + 'ROMEO.01_1d_uo_20181[0-1]*_U.nc'))
vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_201807*_V.nc') + glob(data_path + 'ROMEO.01_1d_vo_20181[0-1]*_V.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

variables = {'U': 'uo',
             'V': 'vo'}

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

pset = ParticleSet.from_line(fieldset=fieldset, 
                             size=200, 
                             pclass=JITParticle,
                             start=(-47.4075912, 1.7378593), 
                             finish=(-47.2618878, 1.5850544), 
                             repeatdt=timedelta(days=1),
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/dmanral/Atlanteco_Romeo_July-Nov2018_OFF250km_z0m.nc", outputdt=timedelta(days=1))

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=300,                       
             output_file=output_file)

output_file.close()
