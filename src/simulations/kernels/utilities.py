from parcels.tools import StatusCode 

def PreventThroughSurfaceError(particle, fieldset, time):
    if particle.depth < fieldset.Surf_Z0:
        particle.depth = fieldset.Surf_Z0


def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        print("Particle [%d] deleted: (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
        particle.delete()

def sudo_AdvectionRK4_3D(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.

    Function needs to be converted to Kernel object before execution.
    """
    (u1, v1, w1) = fieldset.UVW[particle]
    w1 = -w1
    lon1 = particle.lon + u1*.5*particle.dt
    lat1 = particle.lat + v1*.5*particle.dt
    dep1 = particle.depth + w1*.5*particle.dt
    (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1, particle]
    w2 = -w2
    lon2 = particle.lon + u2*.5*particle.dt
    lat2 = particle.lat + v2*.5*particle.dt
    dep2 = particle.depth + w2*.5*particle.dt
    (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2, particle]
    w3 = -w3
    lon3 = particle.lon + u3*particle.dt
    lat3 = particle.lat + v3*particle.dt
    dep3 = particle.depth + w3*particle.dt
    (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3, particle]
    w4 = -w4
    particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
    particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
    particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt


def ParticleSinking(particle, fieldset, time):
    particle.depth += fieldset.sinking_speed * particle.dt # depth adjusted by the speed only.
    