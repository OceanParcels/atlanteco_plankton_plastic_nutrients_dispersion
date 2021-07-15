def SampleTSFields(particle, fieldset, time):
    
    temp = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    sal = fieldset.S[time, particle.depth, particle.lat, particle.lon]
        
    if temp < particle.min_temp:
        particle.min_temp = temp
    elif temp > particle.max_temp:
        particle.max_temp = temp
    
    if sal < particle.min_sal:
        particle.min_sal = sal
    elif sal > particle.max_sal:
        particle.max_sal = sal
    

