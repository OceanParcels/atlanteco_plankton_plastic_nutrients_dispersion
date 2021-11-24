import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode


data_path = '/data/oceanparcels/input_data/NEMO-MEDUSA/ORCA025-N006/'
mesh_mask = data_path + 'domain/coordinates.nc'

ufiles =  sorted(glob(data_path + 'means/ORCA025-N06_2015*d05U.nc'))
vfiles =  sorted(glob(data_path + 'means/ORCA025-N06_2015*d05V.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': vfiles}
            }
variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}}

# filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': ufiles},
#              'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': vfiles}}

# filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
#              'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

# variables = {'U': 'uo',
#              'V': 'vo'}

# # dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}

# dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'},
#               'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}}

simulation_start = datetime(2015, 1, 3, 12, 0, 0)
simulation_end = datetime(2015, 12, 29, 12, 0, 0)

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize=False)

# chs = {'time': ('time_counter', 31), 'depth': ('depthu', 50), 'lat': ('gphif', 3896), 'lon': ('glamf', 1903)}

# fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize=chs)
print(fieldset.U.grid.__dict__)

# stations_data=pd.read_csv(r'UvAStations_Children.csv')
stations_data=pd.read_csv(r'Stations.csv')


pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=stations_data['Longitude'],
                             lat=stations_data['Latitude'],
#                              depth=[0 for i in range(len(stations_data))],
                             time=simulation_start,
                            repeatdt=timedelta(days=5))
                            
output_file = pset.ParticleFile(name="/scratch/dmanral/Orca2015_UvA_StationsRelease_1Yr_rtdt5d_z0m.nc", outputdt=timedelta(days=5))


def delete_particle(particle, fieldset, time):
    particle.delete()
    

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=600,                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()

