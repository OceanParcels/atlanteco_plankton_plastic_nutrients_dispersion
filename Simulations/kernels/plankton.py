import math
# zooplankton_drift kernel: 
# defines the vertical movement of a basic plankton particle that performs Diel Vertical Migration
def ZooplanktonDrift(particle, fieldset, time):    
    
    migration_speed = 0.035     #(m/s)
 
    min_depth = 10      #(m)
    max_depth = 350     #(m) Check this

    # maximum displacement in one dt=0.08 * 300= 24 m (in 5 minutes)
    # https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.1002/lno.10219   (FIGURE 6)
   
    max_displacement = migration_speed * particle.dt
    
    #TODO error handling
    
    # Extract the current hour of the day from the origin timestamp from velocity files(data available from) and time since the origin time.
    
    total_seconds = fieldset.start_time + time
    total_hour = total_seconds/3600
    current_hour = math.fmod(total_hour, 24)
#     current_hour=particle.tod/3600
    
    sunrise_utc = fieldset.Sunrise[time, particle.depth, particle.lat, particle.lon]
    sunset_utc = fieldset.Sunset[time, particle.depth, particle.lat, particle.lon]

#     if Scipy particle, use datetime directly
#     current_datetime = fieldset.start_time + timedelta(seconds=time)
#     current_hour=current_datetime.hour+current_datetime.minute/60 + current_datetime.seconds/3600
        
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
#     print("particleid: %d, time: %f, current_hour: %f,` lat: %f, lon: %f, sunrise: %f, sunset: %f, depth: %f " % (particle.id, time,current_hour, particle.lat, particle.lon, sunrise_utc, sunset_utc, particle.depth))