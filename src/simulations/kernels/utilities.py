import math

def DayoftheYear(particle, fieldset, time):
    #     print("INITI time:%f, tod:%f, doy:%f " % (time, particle.tod, particle.doy))
    # update time of the day
    if particle.tod + particle.dt >= 86400:  # 24 hours completed
        particle.tod = particle.tod - 86400
        # update year if needed and day of the year
        if (particle.is_leap_year and particle.doy + 1 > 366) or (
                not particle.is_leap_year == 0 and particle.doy + 1 > 365):
            particle.doy = 1
            particle.current_year += 1
            # update leap year status
            particle.is_leap_year = math.fmod(particle.current_year, 4) == 0 and (
                    math.fmod(particle.current_year, 100) != 0 or math.fmod(particle.current_year, 400) == 0)
        else:
            particle.doy += 1

    particle.tod += particle.dt
    print("FINAL time:%f, tod:%f, doy:%f " % (time, particle.tod, particle.doy))

def PreventThroughSurfaceError(particle, fieldset, time):
    if particle.depth < 0.5:
        particle.depth = 0.5

def delete_particle(particle, fieldset, time):
    print("Particle [%d] deleted: (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()
