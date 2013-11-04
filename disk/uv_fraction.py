"""
Created on Fri Oct 25 11:39:32 2013

Author: Sam Connolly

uvfraction.py

Calculates the fractional contribution of each radius to the total UV,
and gives the radius within which a given fraction of the UV is emitted
(of the range of radii given) - default UV range of 10-400 nm

"""

from disk import *

def uv_fraction ( m, mdot, rmin, rmax, nfreq = 1000, nrings = 100, fraction = 0.9,
                                w1 = 10, w2 = 400):
        '''
        Function:
            
            uv_fraction creates arrays of radii and, fractional and cumulative 
            contribution to total uv emission (), and gives a radius within which 
            a given fraction of this emission is contained.
        
        Arguments:
            
            m:       
                mass of cental object in msol
            mdot:             
                accretion rate in msol / yr
            rmin, rmax:           
                minimum and maximum radius in cm
            nfreq:             
                number of frequency points [optional]
            nrings:            
                number of disk annuli [optional]
            fraction         : 
                fraction of uv emission to find radius [optional]
            w1, w2:           
                wavelength limits in nm [optional]
                
        Also requires:
                constants.py, functions teff and tdisk from disk.py
                                                                                                                                
        Returns:
                     uv:        
                         (array)    Fractional UV emission at each radius
                     uvCum:     
                         (array)    Cumulative UV emission fraction at each radius
                     rdisk:     
                         (array)    Radii corresponding to these fractions        
                     rFraction: 
                         (float)    Radius within which the given fraction is contained
                                                                
        '''                                        

        # reference temperature of the disk
        tref=tdisk(m, mdot, rmin)
                                                      
        mdot = mdot * MSOL; m = m * MSOL
        
        # change wavelength to frequency, create array of number of 
        # frequencies specified as optional arguments, linear spaced array
        f1,f2 = (C/100.)/(w2*1e-9),(C/100.)/(w1*1e-9)
        freq=np.linspace( f1, f2, nfreq)
                
        # logarithmically spaced radii
        rtemp = np.logspace(np.log10(rmin), np.log10(rmax), num = nrings)

        # uv contribution at each radius, total uv emission                                                        
        uv    = np.zeros(len(rtemp)-1)
        uvTot = 0

        # mean radius of each annulus
        rdisk = np.array([])
       
        # -- loop over annuli ---------------------------------------------------
                                                                
        for j in range(len(rtemp)-1):
                
                # rdisk contains midpoint values for each annulus
                rdisk = np.append(rdisk,(rtemp[j]+rtemp[j+1])/2.0)
                
                # divide by min radius
                r =rdisk[j]/rmin
                
                # area of annulus
                area = PI * (rtemp[j+1]*rtemp[j+1] - rtemp[j]*rtemp[j])
                
                t = ( teff(tref,r) )                # effective temperature of annulus
                                                                                                                        
                for i in range(len(freq)):
                                                                                                                                        
                                      # calculate emission at at each frequency from annulus
                        value = planck_nu(t,freq[i]) * area * PI * 2. 

                        uv[j] += value        # Add to UV emission at this radius
                        uvTot += value   # Add to total UV emission
                                                
        uv /= uvTot # divide UV emission, to give a fraction of total UV
        
        # calculate peak radius for UV emission
                
        rMax  = 0
        uvMax = 0  
        
        for i in range(len(uv)):
            
            if uv[i] > uvMax:
                
                uvMax = uv[i]
                rMax  = rdisk[i]
                     
        print "Radius of peak UV emission is {0} cm = {1} Rg"\
                    .format(rMax,rMax/Schwarz(m/MSOL))             
                     
        #--- Calculate the cumulative UV fraction at each radius ----------------
                                
        uvCum = np.zeros(len(rtemp)-1)                                
        uvCum[0] = uv[0]
        foundFraction = False
                                
        for i in range(1,len(uv)):
                                        
                uvCum[i] = uv[i]+uvCum[i-1] # cumulatively add UV emission fraction
                                
                           # check for desired fraction being reached        
                if uvCum[i] > fraction and foundFraction == False:  
                        foundFraction = True # stop looking for this fraction!
                        rFraction  = rtemp[i]        
                        rFractionS = rFraction /Schwarz(m/MSOL)                                                                                                                
                                          # print result 
                        print "{0}% of UV emission from within {1} cm = {2} Rg"\
                                   .format(fraction*100.,rFraction,rFractionS)
        
        rdisk /= Schwarz(m/MSOL) # change to Rg
                                                                                                                
        # return cumulative uv fraction and radius arrays, fraction radius                
        return uv, uvCum,rdisk,rFraction
