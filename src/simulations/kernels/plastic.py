import math
from parcels import JITParticle, Variable
import numpy as np

'''
Adopted from from Reint and Delphine's work:
https://github.com/OceanParcels/biofouling_3dtransport_2/blob/main/Simulation/kernels.py

'''

class Plastic_Particle(JITParticle):
    '''
    Class declaring properties of a plastic particle and seawater properties needed for the kernel. Assuming a spherical microplastic
    '''
    # Seawater properties
    ambient_density = Variable('ambient_density', dtype=np.float32, to_write=False)    
    
    # Plastic properties
    bouyant_vertical_velocity = Variable('bouyant_vertical_velocity', dtype=np.float32, to_write=True)  # vertical rising or settling velocity due to bouyancy of the particle
    plastic_length = Variable('plastic_length', dtype=np.float32, to_write='once')
    plastic_density = Variable('plastic_density', dtype=np.float32, to_write='once')


# def RK4_Kooi_buoyant_plastic(particle,fieldset,time):  
    # """Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.

    # Function needs to be converted to Kernel object before execution.
    # """
    # (u1, v1, w1) = fieldset.UVW[particle]
    # w1 = -w1
    # lon1 = particle.lon + u1*.5*particle.dt
    # lat1 = particle.lat + v1*.5*particle.dt
    # dep1 = particle.depth + w1*.5*particle.dt
    # (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1, particle]
    # w2 = -w2
    # lon2 = particle.lon + u2*.5*particle.dt
    # lat2 = particle.lat + v2*.5*particle.dt
    # dep2 = particle.depth + w2*.5*particle.dt
    # (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2, particle]
    # w3 = -w3
    # lon3 = particle.lon + u3*particle.dt
    # lat3 = particle.lat + v3*particle.dt
    # dep3 = particle.depth + w3*particle.dt
    # (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3, particle]
    # w4 = -w4
    # particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
    # particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
    # RK4_w = (w1 + 2*w2 + 2*w3 + w4)
    # RK4_dz = (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt

    # """
    # Kernel to compute the vertical velocity (vs) of particles due to their different sizes and densities (no biofouling effect included)- based on Kooi et al. 2017 model 
    # """
    
    # #------ Profiles from MEDUSA or Kooi theoretical profiles -----
    # kin_visc = particle.ambient_kinematic_viscosity # kinematic viscosity[m2 s-1]
    # rho_sw = particle.ambient_density    # seawater density[kg m-3]       

    # #------ Constants -----
    # g = fieldset.gravitational_acc  # gravitational acceleration [m s-2]
    # rho_tot =  particle.plastic_density #  [kg m-3]
    # dn = 2. * (particle.plastic_length)  # equivalent spherical diameter [m]

    # delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    # dstar = ((rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * kin_visc ** 2.)  # [-]
            
    # if dstar > 5e9:
    #     w = 1000.
    # elif dstar < 0.05:
    #     w = (dstar ** 2.) * 1.71E-4
    # else:
    #     w = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
    #                 0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))
     
    # #------ Settling of particle -----
    # if delta_rho > 0:  # sinks
    #     vs = (g * kin_visc * w * delta_rho) ** (1. / 3.)
    # else:  # rises
    #     a_del_rho = delta_rho * -1.
    #     vs = -1. * (g * kin_visc * w * a_del_rho) ** (1. / 3.)  # m s-1

    # # total vertical displacement due to ocean velocities and bouyancy related  
    # total_dz = RK4_dz + vs * particle.dt
    # new_depth = particle.depth + total_dz
    # if new_depth < fieldset.Surf_Z0:  # NEMO's 'surface depth'
    #     # vs = 0
    #     particle.depth = fieldset.Surf_Z0
    # else:
    #     particle.depth = new_depth

    # particle.bouyant_vertical_velocity = vs  # vertical velocity[m s-1] 
    # # print("Particle [%d] bouyant: (%g %g %g %g %g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.ambient_density, particle.bouyant_vertical_velocity, particle.time, RK4_dz, RK4_w))


def RK4_Kooi_vertical_displacement(particle, fieldset, time):
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
    RK4_w = (w1 + 2*w2 + 2*w3 + w4)
    RK4_dz = (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt

    # from Mikael Kaandorp's Code
    # https://github.com/OceanParcels/Global_Analysis_Mikael/blob/main/kernels_v3.py#L128
    #based on Kooi et al.2016
    """
    Calculate settling velocity based on plastic properties (plastic_length, plastic_density),
    and seawater properties (rho_sw, sw_(kin)_visc )
    Assumes spherical particle
    """

    g = fieldset.gravitational_acc  # gravitational acceleration [m s-2]
    temp = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    rho_sw = particle.ambient_density  # [kg m-3]
    
    mu_w = 4.2844E-5 + (1 / ((0.157 * (temp + 64.993) ** 2) - 91.296))
    A = 1.541 + 1.998E-2 * temp - 9.52E-5 * temp ** 2
    B = 7.974 - 7.561E-2 * temp + 4.724E-4 * temp ** 2
    S_sw = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
    ambient_dynamic_viscosity = mu_w * (1 + A * S_sw + B * S_sw ** 2)
    ambient_kinematic_viscosity = ambient_dynamic_viscosity / particle.ambient_density

    r_pl = 0.5 * particle.plastic_length

    # ------ Volumes -----
    rho_tot = particle.plastic_density # total density [kg m-3]

    dn = 2. * (r_pl)  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
    delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    dstar = (math.fabs(rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * ambient_kinematic_viscosity ** 2.)  # [-]

    if dstar > 5e9:
        w_star = 265000
    elif dstar < 0.05:
        w_star = (dstar ** 2.) * 1.71E-4
    else:
        w_star = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                    0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))

    # ------ Settling velocity of particle -----
    if delta_rho > 0:  # sinks
        vs_new = (g * ambient_kinematic_viscosity * w_star * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs_new = -1. * (g * ambient_kinematic_viscosity * w_star * a_del_rho) ** (1. / 3.)  # m s-1
    
    # total vertical displacement due to ocean velocities and bouyancy related  
    total_dz = RK4_dz + vs_new * particle.dt
    new_depth = particle.depth + total_dz
    if new_depth < fieldset.Surf_Z0:  # model's 'minimum surface depth (0.39 m), here set to 0.5 m'
        # vs = 0
        particle.depth = fieldset.Surf_Z0
    else:
        particle.depth = new_depth

    particle.bouyant_vertical_velocity = vs_new  # vertical velocity[m s-1] 
    # print("Particle [%d] bouyant: (%g %g %g %g %g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.ambient_density, particle.bouyant_vertical_velocity, particle.time, RK4_dz, RK4_w))

