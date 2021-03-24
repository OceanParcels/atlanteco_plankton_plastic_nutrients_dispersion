import netCDF4 as nc
import numpy as np
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D
from datetime import timedelta, datetime
from kernels.utilities import PreventThroughSurfaceError

data_path='/data/oceanparcels/input_data/NEMO-MEDUSA/ORCA025-N006/'

mesh_mask = data_path + 'domain/coordinates.nc'

ufiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05U.nc'))
vfiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05V.nc'))
wfiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05W.nc'))

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

simulation_start = datetime(2013, 1, 2, 12, 0, 0)
simulation_end = datetime(2013, 12, 28, 12, 0, 0)

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

time_step = 300

pset = ParticleSet.from_line(fieldset=fieldset, 
                             size=150, 
                             pclass=JITParticle,
                             start=(-47.07351,  1.50464), 
                             finish=(-42.23952, 1.50464), 
                             time=simulation_start, 
                             repeatdt= timedelta(days=1),
                             depth=50) 

output_file_path="/scratch/manra003/Plankton_withoutDVM_3D_1Y_150Prtdt1yr_z50.nc"
output_file = pset.ParticleFile(name=output_file_path, outputdt=timedelta(hours=1))

pset.execute(AdvectionRK4_3D + pset.Kernel(PreventThroughSurfaceError),   
             endtime=simulation_end,             
             dt=time_step,                       
             output_file=output_file)

# pset.repeatdt=None

# pset.execute(AdvectionRK4,   
#              endtime=datetime(2015,12,29,12,0,0),             
#              dt=time_step,                       
#              output_file=output_file)

output_file.close()
      