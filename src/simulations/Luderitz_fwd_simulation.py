import netCDF4 as nc
from parcels import Field, FieldSet, ParticleSet, JITParticle, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.statuscodes import ErrorCode
from datetime import timedelta, datetime, date
import numpy as np
import sys
import parcels
import kernels.utilities as util
import pandas as pd
import kernels.plankton as plankton
import kernels.plastic as plastic
import utils.sunrise_sunset as sun
import dask

dask.config.set({"array.slicing.split_large_chunks": False})
print(parcels.__version__)

if __name__ == '__main__':
    ### ------------------------PRESET VARIABLES------------------------ ###
    dt = timedelta(minutes=10) # integration dt
    output_days= timedelta(days=1) # output dt in days
    start_day = 1
    simulation_dt = 100

    # sinking speed from Pitcher et al. 1989
    # mean Phytoplankton Carbon (PPC) sinking rate for all taxonomic compositions and depths = 0.25 m/d
    # 
    sinking_speed_mpd =  0.25 #

    ### ------------------------PASSED ARGUMENTS------------------------ ###
    # region arguments
    args = sys.argv
    assert len(args) == 6
    start_year = np.int32(args[1]) 
    start_mon = np.int32(args[2])
    mon_name = date(1900, start_mon, 1).strftime('%b')
    release_depth = np.float32(args[3])
    rk_mode = args[4] # types= '3D_BP' '2D' '3D' 'DVM' 'sinking'
    resolution = args[5]

    if rk_mode == '3D_BP':
        plastic_length = 1e-3 # in m 
        plastic_density = 1025   # in kg/m^3
        print("plastic_length :{0} m, plastic_density: {1} kg/m^3".format(plastic_length, plastic_density))

    if rk_mode != '2D':
        min_ind, max_ind = 0, 50 # load all depths if not 2D
    else: # load only specific levels - only mentioned for <1 and 100m 
        if release_depth <= 1:
            min_ind, max_ind = 0, 3
        else:
            raise ValueError('Depth indices have not been setup for 2D release.')
    # endregion

    data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
    project_data_path = '/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'
    mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

    # Model output is given at 12:00:00 H, we will start the simulation at 00:00:00 H
    simulation_start = datetime(start_year, start_mon, start_day, 0, 0, 0)  
    # since we want the release to start at the same depth = 1 m, i.e., plankton is at the surface start at night, one day before the start is also added
    days= [simulation_start-timedelta(days=1)]+[simulation_start+timedelta(days=i) for i in range(simulation_dt+1)]

    ufiles = [data_path + 'ROMEO.01_1d_uo_{0}{1}{2}_grid_U.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
    vfiles = [data_path + 'ROMEO.01_1d_vo_{0}{1}{2}_grid_V.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
    wfiles = [data_path + 'ROMEO.01_1d_wo_{0}{1}{2}.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]

    def define_2Dvariables():
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                    'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles}}

        variables = {'U': 'uo',
                    'V': 'vo'}

        dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
        return filenames, variables, dimensions


    def define_3Dvariables():
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                    'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                    'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}

        variables = {'U': 'uo',
                    'V': 'vo',
                    'W': 'wo'}

        dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
        return filenames, variables, dimensions

    if rk_mode == '2D':
        print('load 2D fields')
        filenames, variables, dimensions = define_2Dvariables()
    else:
        print('load 3D fields')
        filenames, variables, dimensions = define_3Dvariables()

    #'lon': range(605,1903)= 60W,21E  'lat': range(0,1877) = 78S- 0Eq, (0,1877)= 78S to Eq 'depth': range(min_ind, max_ind), 
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': range(min_ind, max_ind), 'lon': range(0,1903), 'lat': range(0,1877)}, chunksize=False) 

    # reversing signs for W, since depth increases 
    # the following doesnt work- for now- changing signs of w in the kernels
    # if rk_mode!='2D':
    #     fieldset.W.set_scaling_factor(-1.0)

    # minimum depth a particle can attain in the simulations.
    fieldset.add_constant('Surf_Z0', 0.5) # in m 

    if rk_mode == 'DVM':
        print("Loading DVM related fieldset information")
        fieldset.add_constant('Plankton_speed', 0.035) # in m/s
        fieldset.add_constant('Plankton_min_depth', 0.5) # in m
        fieldset.add_constant('Plankton_max_depth', 150) # in m

        # store the t0 time- time between the timestamp of first file and the simulation start time.  this is to compute the current hour of any given time during advection.
        # this is a hack 
        time_zero_totalseconds = (datetime.strptime(str(fieldset.U.grid.time_origin)[:19],'%Y-%m-%dT%H:%M:%S') - simulation_start).total_seconds()
        fieldset.add_constant('start_time', time_zero_totalseconds)

        sunrise_da, sunset_da = sun.load_sunrise_sunset(days, project_data_path)
        
        fieldset.add_field(Field.from_xarray(sunrise_da,
                                            'Sunrise', 
                                            dimensions={'lon':'lon','lat':'lat', 'time':'time'},
                                            mesh='spherical'))

        fieldset.add_field(Field.from_xarray(sunset_da,
                                            'Sunset', 
                                            dimensions={'lon':'lon','lat':'lat', 'time':'time'},
                                            mesh='spherical'))

    if rk_mode == '3D_BP':
        print("Loading Bouyant plastic related fieldset information. Radius:{0} m and Density:{1} kg/m^3".format(plastic_length, plastic_density))
        fieldset.add_constant('gravitational_acc', 9.81) # in m/s^2
        # load temperature and salinity fields
        tfiles = [data_path + 'ROMEO.01_1d_thetao_{0}{1}{2}_grid_T.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
        sfiles = [data_path + 'ROMEO.01_1d_so_{0}{1}{2}_grid_T.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
        filenames_BP = {'cons_temperature': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': tfiles},
                        'abs_salinity': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': sfiles}}

        variables_BP = {'cons_temperature': 'thetao',
                        'abs_salinity': 'so'}
        
        dimensions_BP = {'cons_temperature': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'},
                        'abs_salinity': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'}}
        
        BP_fieldset = FieldSet.from_nemo(filenames_BP, variables_BP, dimensions_BP, indices={'depth': range(min_ind, max_ind), 'lon': range(0,1903), 'lat': range(0,1877)}, chunksize=False)

        # For simplification purpose, loading potential temperature from ocean data as conservative temperature(CT)- difference primarily occurs below 1 km depth.
        # ideally convert potential temperature to conservative temperature
        # CT and S needed to compute seawater density. 
        fieldset.add_field(BP_fieldset.cons_temperature)  # degree_C
        fieldset.add_field(BP_fieldset.abs_salinity)      # psu

    if rk_mode == 'sinking':
        fieldset.add_constant('sinking_speed', sinking_speed_mpd/(24*3600)) #speed in m/day to m/s

    u_file = nc.Dataset(ufiles[0])
    ticks = u_file['time_counter'][:][0]
    modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

    assert simulation_start >= modeldata_start

    coords = pd.read_csv(project_data_path + 'Release_points_{0}grid.csv'.format(resolution))

    if release_depth == 0:
        depth_arg = [release_depth for i in range(len(coords['Longitude']))]
    else:
        depth_arg = [release_depth for i in range(len(coords['Longitude']))]


    if rk_mode == '3D_BP':
        # list of radii and densitites for each particle.
        n_length = [plastic_length for i in range(len(coords['Longitude']))]
        n_density = [plastic_density for i in range(len(coords['Longitude']))]
        
        pset = ParticleSet.from_list(fieldset=fieldset, 
                                    pclass=plastic.Plastic_Particle,
                                    lon=coords['Longitude'],
                                    lat=coords['Latitude'],
                                    depth=depth_arg,
                                    time=simulation_start,
                                    plastic_length=n_length,
                                    plastic_density=n_density)
    else:
        pset = ParticleSet.from_list(fieldset=fieldset, 
                                    pclass=JITParticle,
                                    lon=coords['Longitude'],
                                    lat=coords['Latitude'],
                                    depth=depth_arg,
                                    time=simulation_start)
    pset.populate_indices()                            
    output_file = pset.ParticleFile(name="/nethome/manra003/analysis/dispersion/simulations/{0}/BenguelaUpwR_{1}res_{2}{3}_{4}z_{5}days.zarr".format(rk_mode, resolution, mon_name, start_year, int(release_depth), simulation_dt),                               
                                    outputdt=output_days)

    if rk_mode == '3D':
        kernels = pset.Kernel(util.sudo_AdvectionRK4_3D) + pset.Kernel(util.PreventThroughSurfaceError)
    elif rk_mode == 'DVM':  # as of now onyl using uv at z (RK4) and vertical position is dermined buy plankton drift only. - to confirm this!
        kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(plankton.ZooplanktonDrift) 
    elif rk_mode == '3D_BP':
        kernels = pset.Kernel(util.PolyTEOS10_bsq) + pset.Kernel(plastic.RK4_Kooi_vertical_displacement)
    elif rk_mode == 'sinking':
        kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(util.ParticleSinking)
    else:
        kernels= pset.Kernel(AdvectionRK4)

    output_file.add_metadata("Parcels_version", str(parcels.__version__))
    output_file.add_metadata("Days", str(simulation_dt))
    output_file.add_metadata("integration dt in seconds", str(dt))
    output_file.add_metadata("output_dt", str(output_days))
    output_file.add_metadata("time_execution", datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))

    pset.execute(kernels,                
                runtime=timedelta(days=simulation_dt),
                dt= dt,                       
                output_file=output_file,
                recovery={ErrorCode.ErrorOutOfBounds: util.delete_particle})

    output_file.close()