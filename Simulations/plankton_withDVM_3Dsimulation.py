import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D
from datetime import timedelta, datetime
from kernels.plankton import ZooplanktonDrift

data_path='/data/oceanparcels/input_data/NEMO-MEDUSA/ORCA025-N006/'

mesh_mask = data_path + 'domain/coordinates.nc'

ufiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05U.nc'))
vfiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05V.nc'))
wfiles = sorted(glob(data_path + 'means/ORCA025-N06_201501*d05W.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
             'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}
            }

variables = {'U': 'uo',
             'V': 'vo',
             'W': 'wo'}

dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
             }

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')


# extract the origin timestamp from one of the velocity files
u_temp = nc.Dataset(ufiles[0])
ticks = u_temp['time_counter'][:][0]
time_zero = datetime(1900, 1, 1) + timedelta(seconds=ticks)
time_zero_totalseconds = time_zero.hour * 60 * 60 + time_zero.minute * 60 + time_zero.second

# to make this value available to the custom kernel during execution add it to fieldset
fieldset.add_constant('start_time',time_zero_totalseconds)

pset = ParticleSet.from_line(fieldset=fieldset, 
                             size=10, 
                             pclass=JITParticle,
                             start=(-47.07351,  1.50464), 
                             finish=(-42.23952, 1.50464), 
                             time=datetime(2015, 1, 3, 12, 0, 0), 
#                              repeatdt= timedelta(days=1),
                             depth=400) # in m; since start time is noon, if it starts at night, depth = 45 m

# two day simulation therefore, outputdt is small, else per day.

output_file_path="/scratch/manra003/Plankton_withDVM_short_run_2D.nc"
output_file = pset.ParticleFile(name=output_file_path, outputdt=300)

kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(ZooplanktonDrift)
pset.execute(kernels,   
             runtime=timedelta(hours=48),
#              endtime=datetime(2015,1,28,12,0,0),             
             dt=300,                       
             output_file=output_file)

output_file.close()
      