import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from datetime import timedelta, datetime
import pandas as pd
from parcels.tools.statuscodes import ErrorCode

data_path = '/data/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_2018010[1-9]*_U.nc'))
vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_2018010[1-9]*_V.nc'))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

simulation_start = datetime(2018, 1, 1, 12, 0, 0)
# simulation_end = simulation_start + timedelta(days=30)
simulation_end = datetime(2018, 1, 2, 12, 0, 0)

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')

coords = pd.read_csv(r'/scratch/manra003/data/Nemo_H3Release_LatLon_Res6.csv')

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitudes'],
                             lat=coords['Latitudes'],
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/scratch/manra003/FullTara_Res6_Jan2018.nc", outputdt=timedelta((simulation_end-simulation_start).days))

def delete_particle(particle, fieldset, time):
    # remove print later
    print("id: %d, lat: %f, lon: %f, time: %f" % (particle.id, particle.lat, particle.lon, time))    
    particle.delete()

pset.execute(AdvectionRK4,                
             endtime=simulation_end,
             dt=timedelta(hours=3),                       
             output_file=output_file,
            recovery= {ErrorCode.ErrorOutOfBounds:delete_particle})

output_file.close()