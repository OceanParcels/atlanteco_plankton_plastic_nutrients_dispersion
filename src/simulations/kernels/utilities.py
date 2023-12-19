from parcels.tools import StatusCode


def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle_ddepth = fieldset.Surf_Z0
        particle.state = StatusCode.Success
    elif particle.state == StatusCode.ErrorOutOfBounds:
        print("Particle [%d] deleted: (%g %g %g %g)" % (
            particle.id, particle.lon, particle.lat, particle.depth, time))
        particle.delete()


def ParticleSinking(particle, fieldset, time):
    # addtional depth adjustment.
    particle_ddepth += fieldset.sinking_speed * particle.dt
