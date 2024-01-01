import netCDF4 as nc
from parcels import Field, FieldSet, ParticleSet, JITParticle, AdvectionRK4, AdvectionRK4_3D
from datetime import timedelta, datetime, date
import numpy as np
import sys
import parcels
import kernels.utilities as util
import pandas as pd
import kernels.plankton as plankton
import utils.sunrise_sunset as sun
import os
import dask

dask.config.set({"array.slicing.split_large_chunks": False})
print(parcels.__version__)

if __name__ == '__main__':
    ### ------------------------PRE-SET VARIABLES------------------------ ###
    dt = timedelta(minutes=10)  # integration dt
    output_days = timedelta(days=1)  # output dt in days
    start_day = 1
    simulation_days = 100
    min_depth = 1e-3
    print("dt, output_dt_days, start_day, total_days, min_depth", dt, output_days, start_day,simulation_days,min_depth)
    # sinking speed from Pitcher et al. 1989
    # mean Phytoplankton Carbon (PPC) sinking rate for all taxonomic compositions and depths = 0.25 m/day
    #
    sinking_speed_mpd = 0.25

    # DVM settings for intemediate plankton particles, min_depth is fixed to release depth:
    plankton_speed = 0.03  # in m/s
    plankton_max_depth = 150  # in m

    ### ------------------------PASSED ARGUMENTS------------------------ ###
    # region arguments
    args = sys.argv
    assert len(args) == 6
    start_year = np.int32(args[1])
    start_mon = np.int32(args[2])
    mon_name = date(1900, start_mon, 1).strftime('%b')
    release_depth = np.float32(args[3])
    # types= '2D' '3D' 'DVM' 'sinking'
    rk_mode = args[4]
    resolution = args[5]

    if rk_mode != '2D':
        min_ind, max_ind = 0, 50  # load all depths if not 2D
    else:  # load only specific levels -
        # w(0,1,2):0. ,   0.7942803,   1.616721 and u():3.952819e-01, 1.201226e+00, 2.041416e+00
        if release_depth <= 1:
            min_ind, max_ind = 0, 3
        # depth range referred to here with w(5,6) : 4.329307 ,   5.3378363 and u(5,6):4.825891e+00, 5.866306e+00
        elif release_depth == 5:
            min_ind, max_ind = 5, 7
        else:
            raise ValueError(
                'Depth indices have not been setup for 2D release.')
    # endregion

    data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
    project_data_path = '/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'
    mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

    # Model output is given at 12:00:00 H and we want to start the simulation at 00:00:00 H
    simulation_start = datetime(start_year, start_mon, start_day, 0, 0, 0)
    # since we want the release to start at the same depth = 1 m, i.e., plankton is at the surface start at night, one day before the start is also added for interpolation
    days = [simulation_start-timedelta(days=1)]+[simulation_start+timedelta(days=i)
                                                 for i in range(simulation_days+1)]

    ufiles = [data_path + 'ROMEO.01_1d_uo_{0}{1}{2}_grid_U.nc'.format(
        d.strftime("%Y"), d.strftime("%m"), d.strftime("%d")) for d in days]
    vfiles = [data_path + 'ROMEO.01_1d_vo_{0}{1}{2}_grid_V.nc'.format(
        d.strftime("%Y"), d.strftime("%m"), d.strftime("%d")) for d in days]
    wfiles = [data_path + 'ROMEO.01_1d_wo_{0}{1}{2}.nc'.format(
        d.strftime("%Y"), d.strftime("%m"), d.strftime("%d")) for d in days]

    def define_2Dvariables():
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                     'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles}}

        variables = {'U': 'uo',
                     'V': 'vo'}

        dimensions = {'lon': 'glamf', 'lat': 'gphif',
                      'depth': 'depthw', 'time': 'time_counter'}
        return filenames, variables, dimensions

    def define_3Dvariables():
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                     'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                     'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}

        variables = {'U': 'uo',
                     'V': 'vo',
                     'W': 'wo'}

        dimensions = {'lon': 'glamf', 'lat': 'gphif',
                      'depth': 'depthw', 'time': 'time_counter'}
        return filenames, variables, dimensions

    if rk_mode == '2D':
        print('load 2D fields')
        filenames, variables, dimensions = define_2Dvariables()
    else:
        print('load 3D fields')
        filenames, variables, dimensions = define_3Dvariables()

    # 'lon': range(0,1903)= 98W,21E  'lat': range(0,1877) = 78S- 0Eq, 'depth': range(min_ind, max_ind),
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': range(
        min_ind, max_ind), 'lon': range(0, 1903), 'lat': range(0, 1877)}, chunksize=False)

    # minimum depth a particle can attain in the simulations.
    fieldset.add_constant('Surf_Z0', min_depth)  # in m

    if rk_mode == 'DVM':
        print("Loading DVM related fieldset information")
        fieldset.add_constant('Plankton_speed', plankton_speed)  # in m/s
        fieldset.add_constant('Plankton_min_depth', release_depth)  # in m
        fieldset.add_constant('Plankton_max_depth', plankton_max_depth)  # in m

        # store the t0 time- time between the timestamp of first file and the simulation start time.  this is to compute the current hour of any given time during advection.
        # this is a hack
        time_zero_totalseconds = (datetime.strptime(str(fieldset.U.grid.time_origin)[
                                  :19], '%Y-%m-%dT%H:%M:%S') - simulation_start).total_seconds()
        fieldset.add_constant('start_time', time_zero_totalseconds)

        sunrise_da, sunset_da = sun.load_sunrise_sunset(
            days, project_data_path)

        fieldset.add_field(Field.from_xarray(sunrise_da,
                                             'Sunrise',
                                             dimensions={
                                                 'lon': 'lon', 'lat': 'lat', 'time': 'time'},
                                             mesh='spherical'))

        fieldset.add_field(Field.from_xarray(sunset_da,
                                             'Sunset',
                                             dimensions={
                                                 'lon': 'lon', 'lat': 'lat', 'time': 'time'},
                                             mesh='spherical'))

    if rk_mode == 'sinking':
        # speed in m/day to m/s
        fieldset.add_constant('sinking_speed', sinking_speed_mpd/(24*3600))

    u_file = nc.Dataset(ufiles[0])
    ticks = u_file['time_counter'][:][0]
    modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

    assert simulation_start >= modeldata_start

    coords = pd.read_csv(
        project_data_path + 'Benguela_release_points_{0}grid.csv'.format(resolution))

    depth_arg = [release_depth] * len(coords['Longitude'])

    pset = ParticleSet.from_list(fieldset=fieldset,
                                 pclass=JITParticle,
                                 lon=coords['Longitude'],
                                 lat=coords['Latitude'],
                                 depth=depth_arg,
                                #  lat=-34.125,
                                #  lon=18.375,
                                #  depth=release_depth,
                                 time=simulation_start)
    pset.populate_indices()

    output_folder = "/nethome/manra003/analysis/dispersion/simulations/{0}/{1}/".format(
        rk_mode, start_year)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_file = pset.ParticleFile(name=output_folder + "Benguela_{0}_{1}res_{2}-{3}_{4}z_{5}days.zarr".format(rk_mode, resolution, start_year, str(start_mon).zfill(2), int(release_depth), simulation_days),
                                    outputdt=output_days)

    if rk_mode == '3D':
        kernels = [AdvectionRK4_3D, util.CheckOceanBottom]
    elif rk_mode == 'DVM':
        kernels = [AdvectionRK4_3D, plankton.ZooplanktonDrift]
        output_file.add_metadata(
            "dvm min depth in m", str(release_depth))
        output_file.add_metadata(
            "dvm max depth in m", str(plankton_max_depth))
        output_file.add_metadata(
            "migration speed in m/s", str(plankton_speed))
    elif rk_mode == 'sinking':
        kernels = [AdvectionRK4_3D,
                   util.ParticleSinking, util.CheckOceanBottom]
        output_file.add_metadata(
            "sinking speed in m per day", str(sinking_speed_mpd))
    else:
        kernels = [AdvectionRK4]

    kernels += [util.CheckOutOfBounds]
    print(kernels)

    output_file.add_metadata("Days", str(simulation_days))
    output_file.add_metadata("integration dt in seconds", str(dt))
    output_file.add_metadata("output_dt", str(output_days))
    output_file.add_metadata(
        "time_execution", datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))

    pset.execute(kernels,
                 runtime=timedelta(days=simulation_days),
                 dt=dt,
                 output_file=output_file)
