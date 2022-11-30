# sunRunner1D

Welcome to sunRunner1D 0.00002 (2022-11-29). This is a 1D MHD calculation of a CME pulse propagating (and evolving) through an ambient solar wind, state of which can be specified by the user. The version provided here has been compiled on a Mac (and the executable should run on most of the recent versions of the OS, including Ventura). If you plan to run on a Linux computer or Windows, you will need to recompile the Pluto code, which is in the src directory (see the makefile for more information). 

For a run, you need to provide, at a minimum: 

1.  <run name>, e.g. run001 
2. Number of grid points: ‘low’ - 500, ‘medium’ - 1000, ‘large’ - 2000 
3. Outer boundary (in solar radii) - default is ~257 Rs
4. CME parameters at 30Rs 

a. Velocity pulse, default 2000 km/s (dv)
b. Magnetic Field (Bp) pulse, default 600 nT (db)
c. Density pulse, default is 0.0 cm^-3 (dn)
d. Duration of pulse, default is 10 hours (dur)
e. Profile shape for magnetic field pulse (prof): 
   0 = sin^2 (default) 
   1 = half-sine wave

Additionally, you can specify the conditions of the ambient solar wind at the inner boundary (30 Rs). So, if you are trying to match values observed at 1 AU, you may need to adjust them slightly to make them match better. These parameters are:

v0 - speed
rho0 - density
bp0 - transverse magnetic field
t0 - Temperature
    
Note that, currently, the CME has a sin^2 shape but the magnetic field pulse can be set to have a sin shape with half the period so that is has an alternating sign.

All output files are saved in: runs/<run_name>/output. These include diagnostic files (out.txt), which can be checked during run-time as well as the model output (as dbl binary files). Additionally, three image files:

vars_r0.png
vars_r0.png
tracers.png

Show the variables at the inner boundary, the variables at 1 AU, and the evolution of two tracer particles, which are launched at the leading and trailing edge of the CME pulse. 

For help use: 

./sunRunner1D.py -h

In general, to fully specify all parameters, call sunRunner.py with the following arguments: 

./sunRunner1D.py --run <run_name> --grid <grid_size> --r1 <R1> --t0 <t0> --rho0 <rho0> --v0 <v0> --bp0 <bp0> --dv <del_V> --db <del_Bp> --dn <del_rho> --dur <duration> --prof <profile>

As an example, and the one that is included in the repository (as run001), type:

./sunRunner1D.py --run run001 --grid medium --r1 260 --dv 4000 --db 500 --dn 0 --dur 12 --prof 0
