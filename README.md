# sunRunner1D

Welcome to sunRunner1D. This is a 1D MHD calculation of a CME pulse propagating (and evolving) through an ambient solar wind, the state of which can be specified by the user. The version provided here has been compiled on a Mac (and the executable should run on most of the recent versions of the OS, up to, and including Sonoma 14.2.1). If you plan to run on a Linux computer or Windows, you will need to recompile the Pluto code, which is in the src directory (See instructions below). 

For a run, you need to provide, at a minimum: 

1.  run name, e.g. run001 
2. Number of grid points: ‘low’ - 500, ‘medium’ - 1000, ‘high’ - 2000 
3. Outer boundary (in solar radii) - default is 260 Rs
4. CME parameters at 21.5Rs
   
   a. Velocity pulse, default 500 km/s (dv)

   b. Magnetic Field (Bp) pulse, default 600 nT (db)

   c. Density pulse, default is 0.0 cm^-3 (dn)

   d. Duration of pulse, default is 10 hours (dur)

   e. Profile shape for magnetic field pulse (prof): 

      0 = sin^2 (default) 

      1 = half-sine wave

Additionally, you can specify the conditions of the ambient solar wind at the inner boundary (21.5 Rs). So, if you are trying to match values observed at 1 AU, you may need to adjust them slightly to make them match better. These parameters are:

v0 - speed
rho0 - density
bp0 - transverse magnetic field
t0 - Temperature
    
Note that, currently, the CME has a sin^2 shape but the magnetic field pulse can be set to have a sin shape with half the period so that it has an alternating sign.

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


PLUTO INSTALLATION

The sunRunner1D tool is based on the PLUTO code. 

Hence, PLUTO must first be downloaded and installed from here:

http://plutocode.ph.unito.it/

If you are not familiar with the PLUTO code we recommend that you spend time reading the manual and running the examples. If you would like to understand how the boundary conditions (BC) can be defined please see chapter 5 of the PLUTO manual.  


COMPILING PLUTO FOR sunRunner1D

We have made modifications to only two of the PLUTO routines: init.c and userdef_output.c These routines along with the required definitions.h file are in the 'src' sub-directory.  To compile the PLUTO code for sunRunner1D please go to the src sub-directory:

% cd src

Invoke the PLUTO setup.py script (see Section 1.3 in the PLUTO manual) using:

%  python $PLUTO_DIR/setup.py

Once you are done with this configuration setup you will see a makefile and sysconf.out file that are specific to your computer platform (we recommend selecting the gcc compiler in setup.py) and you can now compile PLUTO:

%make 

Clean the directory:

%make clean

And go back to the directory with the sunRunner1D.py script

%cd ..

You are now ready to use the script.

Please note that the jupyter notebook in the repository is currently not supported and not up-to-date. 
   
For help, please email us at: pete@predsci.com, mbennun@predsci.com
   
