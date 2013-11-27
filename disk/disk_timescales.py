"""
Created on Wed Oct 30 11:47:35 2013

Author: Sam Connolly

Calculate various timescales at a given radius in the disk, calculated using
the thin disc model, as derived by Frank, King and Raine (2002) - "Accretion
Power in Astrophysics"

Requires:
    disk, disky_Const

"""

from disk import *
from disk_const import *
import numpy as np

# ------------------------------------------------------------------------------

def T_light_travel(M,R):

    '''
    Function:
        Calculates the light travel time from a given gravitational radius 
        to the BH (1 Gr)
        
    Args:
        M:
            BH mass in Msolar
        R:
            Disc radius
            
    Returns:
        Light travel time in seconds
    '''

    Rg = Schwarz(M)
    
    R = R * Rg
    
    Rdiff = R-Rg
    
    T = Rdiff/C
    
    return T


# Dynamical timescale

def T_dynamical(M,R):
    
    '''
    Function:
        Calculates dynamical timescale at a given radius of an accretion disk
        
    Args:
            M:           
                BH mass in Msolar
            R:           
                Radius in gravitational radii
                
    Returns:
            Dynamical timescale in seconds
    '''
    
    R10 = Rg/1e10                # radius in 1 x 10 ^10 cm
   
    Tdyn = (100.* M**(-0.5)*R10**(3./2.))

    return Tdyn
    
# Thermal timescale

def T_thermal(M,R,alpha):
     
    '''
    Function:
        Calculates viscous timescale at a given radius of an accretion disk
        
    Args:
            alpha:       
                Viscosity parameter (0-1.0)
            M:           
                BH mass in Msolar
            R:           
                Radius in gravitational radii
                
    Returns:
            Viscous timescale in seconds
    '''
      
    Tth = (T_dynamical(M,R)/alpha)

    return Tth
    
# viscous timescale

def T_viscous(alpha,Mdot,M,R):
    
    '''
    Function:
        Calculates viscous timescale at a given radius of an accretion disk
        
    Args:
            alpha:       
                Viscosity parameter (0-1.0)
            Mdot:        
                Mass accretion rate in Msolar/year
            M:           
                BH mass in Msolar
            R:           
                Radius in gravitational radii
                
    Returns:
            Viscous timescale in seconds
    '''
    
    Mdot16 =  (Mdot*YR*MSOL)/1e16  # mass accretion rate in 1x10^16 g/s
    
    Rg   = R * Schwarz(M)      # radius in cm

    R10 = Rg/1e10                # radius in 1 x 10 ^10 cm
    
    Tvis = (3e5 * (alpha**(-0.8)) * (Mdot16**(-0.3)) * (M**(0.25)) * (R10**(1.25)))

    return Tvis



# Viscous timescale between two points on the disk

def Diff_T_viscous(R1,R2,alpha,Mdot,M):

    '''
    Function:
        Calculates the viscous timescale between two radii of an accretion disk
        
    Args:
            R1:          
                Outer radius in gravitational radii
            R2:          
                Inner radius in gravitational radii
            alpha:       
                Viscosity parameter (0-1.0)
            Mdot:        
                Mass accretion rate in Msolar/year
            M:           
                BH mass in Msolar
                
    Returns:
            Viscous timescale in seconds
    '''
    
    return T_viscous(alpha,Mdot,M,R1) - T_viscous(alpha,Mdot,M,R2)

# Wind Structural change timescale

def T_wind(Rin,v_wind = 1e3):

    '''
    Function:
        Calculates the characteristic timescale for changes in the structure
        of a disc wind with a given innner radius and outflow velocity
        
    Args:
            Rin:          
                Inner wind radius in cm
            v_wind:
                Wind outflow velocity in km s^-1 [optional]
                    
    Returns:
            Wind timescale in seconds
    '''
    
    return ((Rin/1e5)/v_wind)
    
def T_accrete(alpha,Mdot,M,R,nsteps = 10000):

    '''
    Function:
        Calculates the integrated time for matter to be accreted through
        an accretion disc from a given radius, by calculating the inward
        velocity at each radius. Should be about the same as the viscous
        timescale.
        
    Args:
            alpha:          
                The disc viscosity parameter (0-1)
            Mdot:
                Mass accretion rate in Msol per year
            M:
                BH mass in Msolar
            R:
                Radius at which matter starts
            nsteps:
                Number of radii to integrate over [optional]
                    
    Returns:
            Wind timescale in seconds
    '''

    Rg = Schwarz(M)
    Rg10 = Rg / 1e10

    R10 = (R * Rg) / 1e10 # radius in 10^10 cm
    
    Mdot16 = (Mdot*MSOL*YR) / 1e16 # accretion rate in 10^16 g/s

    
    rtemp = np.logspace(np.log10(Rg10), np.log10(R10), num = nsteps)
    
    time = 0

    for i in range(len(rtemp)-1):

        dr = rtemp[i+1]-rtemp[i]
        r10 = (rtemp[i+1]+rtemp[i])/2.0        
        f = (1 - (1./(Rg/r10))) # whatever f is...
        
        Vr = 2.7e4 * alpha**(0.8) * Mdot16**(0.3) * M**(-0.25) * r10**(-0.25) * f**(-14./5.)    
        
        dt = (dr*1e10)/Vr
        
        time += dt
        
    return time
