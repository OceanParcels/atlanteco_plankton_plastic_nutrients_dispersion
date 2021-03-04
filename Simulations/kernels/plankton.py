import math

# defines the vertical movement of a simple plankton particle that performs Diel Vertical Migration(DVM)
def ZooplanktonDrift(particle, fieldset, time):    
    
    migration_speed = 0.035     #(m/s)
 
    min_depth = 10      #(m)
    max_depth = 350     #(m)

    max_displacement = migration_speed * particle.dt
    
    #TODO error handling
    
    # Extract the current hour of the day from the origin timestamp from velocity files(data available from) and time since the origin time.
    
    total_seconds = fieldset.start_time + time
    total_hour = total_seconds/3600
    current_hour = math.fmod(total_hour, 24)
   
    sunrise_utc = fieldset.Sunrise[time, particle.depth, particle.lat, particle.lon]
    sunset_utc = fieldset.Sunset[time, particle.depth, particle.lat, particle.lon]
    
    # Polar locations- Not tested yet:
    #sunrise_utc=-1 or sunset_utc=-2: sun never sets here
    #sunrise_utc=-2 or sunset_utc=-1: sun never rises here
    
    if (current_hour >= sunrise_utc and current_hour < sunset_utc) or (sunrise_utc ==  -1 and sunset_utc ==-2) : #day
        if((particle.depth + max_displacement)>max_depth):
            particle.depth = max_depth
        else:
            particle.depth += max_displacement
            
    else:   #night
        if (particle.depth - max_displacement)<min_depth:
            particle.depth = min_depth
        else:
            particle.depth -= max_displacement
