# Dark-Matter-Subs-and-Globs

## Process:  <br />
out2gyr uses the Python package _clustertools_ to generate six executable files. When executed, the first and second files created by out2gyr converts the ascii initial conditions of, 1 - an nbody star cluster that resides in a smooth galactic potential and 2 - an nbody star cluster that resides in a smooth galactic potential and a moving dark matter substructure potential, to a binary format that is ready to be integrated with the nbody integration code _gyrfalcON_ (GalaxY simulatoR using falcON).  When executed, the third and forth files created by out2gyr begin the _gyrfalcON_ integration of both star cluster environments with the variables passed to the original out2gyr call. After the _gyrfalcon_ integrations are finished fifth and sixth files can be executed to convert back from binary to ascii format to be analyzed. *out2gyr requires functions contained within subhalo_populations.py for the creation of the substructure populations.*

## Background: <br />
This work expands upon Dark Matter Subs and Tracers with the goal of providing a definitive answer to the questions: how much do dark matter subhalos effect the evolution of star clusters, and are star clusters affected in ways in which we can measure.
