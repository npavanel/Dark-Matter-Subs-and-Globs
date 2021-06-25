#imports
from galpy.orbit import Orbit
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.constants as c
import astropy.units as u
from galpy.potential import LogarithmicHaloPotential
from galpy import potential as gp
from galpy.util import bovy_conversion
from galpy.potential import nemo_accname, nemo_accpars
import time as pytime

#GALPY scaling variables
ro=8.
vo=220.
to=bovy_conversion.time_in_Gyr(ro=ro,vo=vo)
mo=bovy_conversion.mass_in_msol(ro=ro,vo=vo)

#gyfalcon units
gyr_m=222288.4543021174 #mass conversion
gyr_v=0.9777922216731284 #velocity conversion
gyr_s=0.0015 #softenting

#integration time
torb = np.linspace(0,12.,250)

#radius and mass at each radius in each potential
r=np.linspace(0,100.,500)

#probability density function
def rndm(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b
    """
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

#define a funciton that returns a population of subhalos 
def subpop_gyr(smf=float, subhalo_mass=None):
    '''
    this funciton returns an array of subhalos that are ready to be passed to gyr falcon for integration.
    
    pop_info holds information about the subhalo population. In order: total number of subhalos, number of 1e6-1e7 subs, number of 1e7-1e8 subs, number of 1e8-1e9 subs
    '''
        
    #initialize arrays for output
    pop_info=[]
    
    #substructure mass fraction
    smf = smf

    #potential creation (total,smooth,clumpy)
    mw = gp.MWPotential2014()
    tp = gp.LogarithmicHaloPotential(amp=1., ro=ro,vo=vo)
    sp = gp.LogarithmicHaloPotential(amp=1.-smf, ro=ro,vo=vo)
    cp = gp.LogarithmicHaloPotential(amp=smf, ro=ro,vo=vo)

    #radius and mass at each radius in each potential
    r=np.linspace(0,100.,5000)
    tmar = []
    smar = []
    cmar = []
    for i in r:
        tmar.append(tp.mass(i/ro))
        smar.append(sp.mass(i/ro))
        cmar.append(cp.mass(i/ro))

    #generate subhalo masses
    if subhalo_mass is None:
        #generate random subhalo masses
        masses = rndm(10**6,10**9,-1,size=1000000)
    else:
        masses=np.ones(int(cmar[-1]/subhalo_mass))*subhalo_mass

    #determine how many subhalos will make up the distribution
    cum_masses=np.cumsum(masses)
    indx=cum_masses<cmar[-1]
    count=np.sum(indx)
    sh_masses = masses[0:count]
    
    #append to output list
    pop_info.append(count)

    #populate the subhalo mass list
    sh_masses = masses[0:count]

    #determine the distribution of subhalos over the mass range
    mass_range=[10**6,10**7,10**8,10**9]
    mass_lower=[10**6,10**7,10**8]
    mass_upper=[10**7,10**8,10**9]
    for i in range(len(mass_lower)):
        indx=(sh_masses >= mass_lower[i])*(sh_masses<mass_upper[i])
        #append to the output list
        pop_info.append(np.sum(indx))

    #generate random subhalo galactocentric radii following the mass profile of the logarithmic halo
    ran=np.random.rand(count)
    rad=np.linspace(r[0],r[-1],count)
    menc=[]
    for i in rad:
        menc.append(cp.mass(i/ro))
    menc=np.array(menc)
    sh_radii=np.interp(ran, menc/menc[-1], rad)

    #initialize the hernquist potentials for the subhalos
    sh_plummer_r=[]
    for i in sh_masses:
        sh_plummer_r.append(1.62*(i/(10**(8)))**(0.5))

    sh_pots=[]
    for i in range(len(sh_masses)):
        sh_pots.append(gp.PlummerPotential(sh_masses[i]/mo,sh_plummer_r[i]/ro,ro=ro,vo=vo))

    #determine the circulr velocity of each subhalo
    sh_circv=[]
    for i in sh_radii:
        sh_circv.append(tp.vcirc(i/ro))  
    sh_circv=np.array(sh_circv)

    #generate random initial conditions for subhalos
    shthea=np.arccos(1.-2.*np.random.rand(len(sh_radii)))
    shphi=2*np.pi*np.random.rand(len(sh_radii))
    shR = sh_radii*np.sin(shthea)
    shZ = sh_radii*np.cos(shthea)

    #switch random ics to cartesian for output file
    shx = shR * np.cos(shphi)
    shy = shR * np.sin(shphi)
    shz = shZ
    shvx = (1/(np.sqrt(3)))*np.random.normal(0,sh_circv)
    shvy = (1/(np.sqrt(3)))*np.random.normal(0,sh_circv)
    shvz = (1/(np.sqrt(3)))*np.random.normal(0,sh_circv)

    #get a column stack for output for nbody
    shorb_nbody=np.column_stack([shx,shy,shz,shvx,shvy,shvz,sh_masses]) 

    return shorb_nbody, sh_plummer_r, pop_info, sp

def cluster_orbit(smf, r):
    
    sp = gp.LogarithmicHaloPotential(amp=1.-smf, ro=ro,vo=vo) # define the smooth potential 
    
    orbit = [r/ro,0,sp.vcirc(r/ro)/vo,0,0,0]
    orbit = Orbit(orbit,ro=ro,vo=vo)
    
    return orbit
    
    
