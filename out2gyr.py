# import: regulars
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from time import time

# clustertools
import clustertools as ctools

# galpys
from galpy.orbit import Orbit
from galpy.potential import LogarithmicHaloPotential
from galpy import potential as gp
from galpy.util import bovy_conversion

# self defined functions that return subhalo population ready for gyr, and an orbit for the cluster
from subhalo_populations import subpop_gyr, cluster_orbit

# GALPY scaling variables
ro=8.
vo=220.
to=bovy_conversion.time_in_Gyr(ro=ro,vo=vo)
mo=bovy_conversion.mass_in_msol(ro=ro,vo=vo)

#integration time
torb = np.linspace(0,12.,1000)

gyr_s=0.00015 # softenting value for stars

# want to define a function that outputs everything
def out2gyr(cluster_rs, num_stars, hmr, smf, kmax, tstop, softenting, step, subhalo_mass=False):
    """
    This function takes a radii for a cluster orbit, a number of stars to include in the cluster, 
    a substructure mass fraction for the potential, and a choice of mass specturm or not for the subhalo population.
    """
    
    # open files for executing on server
    if subhalo_mass is False:
        a2s_file = open('smf'+str(smf)+'_a2s_calls.txt','w')
        gyr_file = open('smf'+str(smf)+'_gyr_calls.txt','w')
        s2a_file = open('smf'+str(smf)+'_s2a_calls.txt','w')
    else:
        a2s_file = open('smf3_'+str(np.log10(subhalo_mass))+'submass_a2s_calls.txt','w')
        gyr_file = open('smf3_'+str(np.log10(subhalo_mass))+'submass_gyr_calls.txt','w')
        s2a_file = open('smf3_'+str(np.log10(subhalo_mass))+'submass_s2a_calls.txt','w')
    
    # loop through a list of radii to produce all clusters at once
    for cluster_r in cluster_rs:
        cluster_orb = cluster_orbit(smf, cluster_r)     # compute a circular orbit for the cluster to orbit at

        # compute a star cluster with num_stars on a circular orbit at cluster_r - base has no subs
        base_cluster = ctools.setup_cluster('limepy', g=1, phi0=5.0, M=float(num_stars), rh=hmr, N=num_stars, orbit=cluster_orb) 
        base_cluster.to_galaxy()
        base_cluster.to_kpckms()

        # compute a star cluster with num_stars on a circular orbit at cluster_r - ptrb has subs
        ptrb_cluster = ctools.setup_cluster('limepy', g=1, phi0=5.0, M=float(num_stars), rh=hmr, N=num_stars, orbit=cluster_orb)
        ptrb_cluster.to_galaxy()
        ptrb_cluster.to_kpckms()

        # produce a mass specturm subhalo population if False, a single mass population if True
        if subhalo_mass is False:
            subs, subs_softenings, sub_info, sp, tp = subpop_gyr(smf, subhalo_mass=None) 
        else:
            subs, subs_softenings, sub_info, sp, tp = subpop_gyr(smf, subhalo_mass) 

        # add the subhalo population to the cluster
        ptrb_cluster.add_stars(subs[:,0],subs[:,1],subs[:,2],subs[:,3],subs[:,4],subs[:,5],subs[:,6])

        # create a list of softening values for every star and append to star softenings
        base_softenings = np.ones(num_stars) * gyr_s             
        ptrb_softenings = np.append(base_softenings, subs_softenings)

        # produce inital condition files in gyrfalcon format, names differing for mass specturm subhalo pop vs single mass pop
        if subhalo_mass is False:
            base_ascii_filename = 'base_IC_' + str(cluster_r) + 'kpc_' + str(int(smf*100)) + 'smf_mspec_' + str(num_stars) + 'stars'
            ptrb_ascii_filename = 'ptrb_IC_' + str(cluster_r) + 'kpc_' + str(int(smf*100)) + 'smf_mspec_' + str(num_stars) + 'stars'
            ctools.util.output.gyrout(base_cluster, filename = base_ascii_filename, eps=base_softenings, epsunits=None, ro=8.0)
            ctools.util.output.gyrout(ptrb_cluster, filename = ptrb_ascii_filename, eps=ptrb_softenings, epsunits=None, ro=8.0)
            #print('Routine Completed')
        else:
            base_ascii_filename = 'base_IC_' + str(cluster_r) + 'kpc_' + str(int(smf*100)) + 'smf_' + str(int(np.log10(subhalo_mass))) + 'submass_' + str(num_stars) + 'stars'
            ptrb_ascii_filename = 'ptrb_IC_' + str(cluster_r) + 'kpc_' + str(int(smf*100)) + 'smf_' + str(int(np.log10(subhalo_mass))) + 'submass_' + str(num_stars) + 'stars'
            ctools.util.output.gyrout(base_cluster, filename = base_ascii_filename, eps=base_softenings, epsunits=None, ro=8.0)
            ctools.util.output.gyrout(ptrb_cluster, filename = ptrb_ascii_filename, eps=ptrb_softenings, epsunits=None, ro=8.0)
            #print('Routine Completed')

        # write to files for execution
        if subhalo_mass is False:
            a2s_file.write('a2s in=' + base_ascii_filename + ' out=BN_' + base_ascii_filename + ' N=' + str(num_stars) + ' read=mxve \n')
            a2s_file.write('a2s in=' + ptrb_ascii_filename + ' out=BN_' + ptrb_ascii_filename + ' N=' + str(sub_info[0] + num_stars) + ' read=mxve \n')
            gyr_file.write('gyrfalcON in=BN_' + base_ascii_filename + ' out=base_EV' + base_ascii_filename[7:] + ' give=mxvpqael tstop=' + str(tstop) + ' eps=-' + str(softenting) + ' step=' + str(step) + ' kmax=' + str(kmax) + ' Nlev=10 fac=0.01 theta=0.6 accname=LogPot accpars=' + str(tp.nemo_accpars(220.,8.)) + ' > gyr_call_' + base_ascii_filename + '.dat 2>&1 & \n')
            gyr_file.write('gyrfalcON in=BN_' + ptrb_ascii_filename + ' out=ptrb_EV' + ptrb_ascii_filename[7:] + ' give=mxvpqael tstop=' + str(tstop) + ' eps=-' + str(softenting) + ' step=' + str(step) + ' kmax=' + str(kmax) + ' Nlev=10 fac=0.01 theta=0.6 accname=LogPot accpars=' + str(sp.nemo_accpars(220.,8.)) + ' > gyr_call_' + ptrb_ascii_filename + '.dat 2>&1 & \n')
            s2a_file.write('s2a in=base_EV' + base_ascii_filename[7:] + ' out=base_EV' + base_ascii_filename[7:] + '.dat \n')
            s2a_file.write('s2a in=ptrb_EV' + ptrb_ascii_filename[7:] + ' out=base_EV' + ptrb_ascii_filename[7:] + '.dat \n')
        else:
            a2s_file.write('a2s in=' + base_ascii_filename + ' out=BN_' + base_ascii_filename + ' N=' + str(num_stars) + ' read=mxve \n')
            a2s_file.write('a2s in=' + ptrb_ascii_filename + ' out=BN_' + ptrb_ascii_filename + ' N=' + str(sub_info[0] + num_stars) + ' read=mxve \n')
            gyr_file.write('gyrfalcON in=BN_' + base_ascii_filename + ' out=base_EV' + base_ascii_filename[7:] + ' give=mxvpqael tstop=' + str(tstop) + ' eps=-' + str(softenting) + ' step=' + str(step) + ' kmax=' + str(kmax) + ' Nlev=10 fac=0.01 theta=0.6 accname=LogPot accpars=' + str(tp.nemo_accpars(220.,8.)) + ' > gyr_call_' + base_ascii_filename + '.dat 2>&1 & \n')
            gyr_file.write('gyrfalcON in=BN_' + ptrb_ascii_filename + ' out=ptrb_EV' + ptrb_ascii_filename[7:] + ' give=mxvpqael tstop=' + str(tstop) + ' eps=-' + str(softenting) + ' step=' + str(step) + ' kmax=' + str(kmax) + ' Nlev=10 fac=0.01 theta=0.6 accname=LogPot accpars=' + str(sp.nemo_accpars(220.,8.)) + ' > gyr_call_' + ptrb_ascii_filename + '.dat 2>&1 & \n')
            s2a_file.write('s2a in=base_EV' + base_ascii_filename[7:] + ' out=base_EV' + base_ascii_filename[7:] + '.dat \n')
            s2a_file.write('s2a in=ptrb_EV' + ptrb_ascii_filename[7:] + ' out=base_EV' + ptrb_ascii_filename[7:] + '.dat \n')
        
    # close files
    a2s_file.close()
    gyr_file.close()
    s2a_file.close()
