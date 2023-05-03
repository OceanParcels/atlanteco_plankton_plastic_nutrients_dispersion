import netCDF4 as nc
import numpy as np
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, Field,Variable
from parcels.tools.converters import TimeConverter
from datetime import timedelta, datetime
from kernels.plankton import ZooplanktonDrift
import math

data_path='/data/oceanparcels/input_data/NEMO-MEDUSA/ORCA025-N006/'

mesh_mask = data_path + 'domain/coordinates.nc'

ufiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05U.nc'))
vfiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05V.nc'))
wfiles = sorted(glob(data_path + 'means/ORCA025-N06_2013*d05W.nc'))

sunrisefile = '/scratch/manra003/SunriseTime_5x5_7d_Atlantic_Jan2015.nc'
sunsetfile = '/scratch/manra003/SunsetTime_5x5_7d_Atlantic_Jan2015.nc'
sunrise_nc = nc.Dataset(sunrisefile,'r')
sunset_nc = nc.Dataset(sunsetfile,'r')

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles}
            }

variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
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

# total number of seconds in that "day" from which the data is available.
time_zero_totalseconds = modeldata_start.hour * 60 * 60 + modeldata_start.minute * 60 + modeldata_start.second
fieldset.add_constant('start_time', time_zero_totalseconds)

time_origin = TimeConverter(np.datetime64(nc.num2date(sunrise_nc['time'][0],'seconds since '+ str(simulation_start.year) + '-01-01',sunrise_nc['time'].calendar)))

fieldset.add_field(Field("Sunrise",
                        data=sunrise_nc['sunrise'][::],
                        lon=sunrise_nc['lon'][:],
                        lat=sunrise_nc['lat'][:],
                        time=sunrise_nc['time'][:],
                        time_origin=time_origin,
                        transpose=False,
#                         time_periodic=
                        interp_method='linear'))

fieldset.add_field(Field("Sunset",
                        data=sunset_nc['sunset'][::],
                        lon=sunset_nc['lon'][:],
                        lat=sunset_nc['lat'][:],
                        time=sunset_nc['time'][:],
                        time_origin=time_origin,
                        transpose=False,
#                         time_periodic=
                        interp_method='linear'))

time_step = 300

pset = ParticleSet.from_line(fieldset=fieldset,  
                             size=150, 
                             pclass=JITParticle,
                             start=(-47.07351,  1.50464), 
                             finish=(-42.23952, 1.50464), 
                             time=simulation_start, 
                             repeatdt= timedelta(days=1),
                             depth=350) 

output_file_path ="/scratch/manra003/Plankton_withDVM_2D_1Y_150Prtdt1yr_z350.nc"
output_file = pset.ParticleFile(name=output_file_path, outputdt=timedelta(hours=1))

kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(ZooplanktonDrift)
pset.execute(kernels,   
             endtime=simulation_end,             
             dt=time_step,                       
             output_file=output_file)

output_file.close()
       