from parcels.tools import StatusCode


def CheckOutOfBounds(particle, fieldset, time):
    # if particle.depth + particle_ddepth > fieldset.Plankton_max_depth:  
    # if particle.state == StatusCode.ErrorThroughSurface:
    #     # print("Particle [%d] moved: (%g %g %g %g %g %g %g)" % (
    #     #     particle.id, particle.lon, particle.lat, particle.depth, particle_dlon, particle_dlat, particle_ddepth, time))
    #     particle_ddepth = 0
    #     particle.state = StatusCode.Success
    #     # particle.depth = fieldset.Surf_Z0

    if particle.depth + particle_ddepth < fieldset.Surf_Z0:
        particle_ddepth = 0
        particle.depth = fieldset.Surf_Z0
    if particle.state == StatusCode.ErrorOutOfBounds:
        print("Particle [%d] deleted: (%g %g %g %g)" % (
            particle.id, particle.lon, particle.lat, particle.depth, time))
        particle.delete()


def ParticleSinking(particle, fieldset, time):
    # additional depth adjustment.
    particle_ddepth += fieldset.sinking_speed * particle.dt


def CheckOceanBottom(particle, fieldset, time):
    # if new location leads to below ocean bottom(previous particle location is below the ocean- positive), particle_dlat and dlon are not updated= 0
    if particle.state < 50 and particle_dlat == 0.0 and particle_dlon == 0.0 and particle.depth > 0:
        particle.delete()
        print("Particle [%d] deleted at bottom: (%g %g %g %g)" % (
            particle.id, particle.lon, particle.lat, particle.depth, time))
