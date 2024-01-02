import math

# defines the vertical movement of a simple plankton particle that performs Diel Vertical Migration(DVM)


def ZooplanktonDrift(particle, fieldset, time):

    # TODO: increasing the speed of plankton from 0 to max as it leaves a given depth
    # and decreasing the speed from max to 0 as the plankton approaches desired depth
    # approximate values for Copepods

    max_vertical_displacement = fieldset.Plankton_speed * particle.dt

    # Extract the current hour of the day from the origin timestamp from velocity files
    # (data available from) and time since the origin time.
    total_seconds = time - fieldset.start_time
    total_hour = total_seconds/3600
    current_hour = math.fmod(total_hour, 24)

    sunrise_utc = fieldset.Sunrise[time,
                                   particle.depth, particle.lat, particle.lon]
    sunset_utc = fieldset.Sunset[time,
                                 particle.depth, particle.lat, particle.lon]

    # Polar locations- Not tested yet:
    # sunrise_utc=-1 or sunset_utc=-2: sun never sets here
    # sunrise_utc=-2 or sunset_utc=-1: sun never rises here

    # additional conditons -maybe for polar: or (sunrise_utc ==  -1 and sunset_utc ==-2) : #day
    if (current_hour >= sunrise_utc and current_hour < sunset_utc):
        # keep going down or stay
        if particle_dlat != 0.0 or particle_dlon != 0.0:
            if ((particle.depth + particle_ddepth + max_vertical_displacement) > fieldset.Plankton_max_depth):
                particle.depth = fieldset.Plankton_max_depth
                particle_ddepth = 0
            else:
                particle_ddepth += max_vertical_displacement

    else:  # night
        if (particle.depth + particle_ddepth - max_vertical_displacement) < fieldset.Plankton_min_depth:
            particle.depth = fieldset.Plankton_min_depth
            particle_ddepth = 0

        else:
            particle_ddepth -= max_vertical_displacement

    # print("FINAL lat:%f, lon:%f, depth:%f, current_hour:%f, sunrise_utc:%f, sunset_utc:%f" % (particle.lat, particle.lon, particle.id, particle.depth, current_hour, sunrise_utc, sunset_utc))
