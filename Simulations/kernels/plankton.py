import math

# zooplankton_drift kernel: 
# defines the vertical movement of a basic plankton particle that performs Diel Vertical Migration
def ZooplanktonDrift(particle, fieldset, time):    
    
    migration_speed = 0.08     #(m/s)
 
    min_depth = 45      #(m)
    max_depth = 400     #(m) Check this

    # maximum displacement in one dt=0.08 * 300= 24 m (in 5 minutes)
    # https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.1002/lno.10219   (FIGURE 6)
   
    max_displacement= migration_speed * particle.dt
    
    #TODO
    # day and night timings are hard coded- works well for equatorial waters
    # maybe use a function on latitude and time to calculate the daylight and night timings 
    # error handling
    
    # day from 6 am to 6 pm 
    # night from 6 pm to 6 am 

    # Extract the current hour of the day from the origin timestamp from velocity files(data available from)
    # and time since the origin time.
    
    total_seconds=fieldset.start_time + time
    total_hour=total_seconds/3600
    current_hour= math.fmod(total_hour, 24)
        
    if current_hour>=6 and current_hour < 18: #day
        if((particle.depth + max_displacement)>max_depth):
            particle.depth=max_depth
        else:
            particle.depth+=max_displacement
            
    else:   #night
        if (particle.depth - max_displacement)<min_depth:
            particle.depth=min_depth
        else:
            particle.depth-=max_displacement
    