# rebound_runs

end_aleksey_config_b.py: Initializes and runs a rebound simulation starting from a config file 
Usage: python end_aleksey_config_b.py --config CONFIG_FILE_NAME (if CONFIG_FILE_NAME is not specified the it default to "config").

The config file is structured as follows

[params]
....

[component 1]
....

[component 2]
.....

The params section is always required and specified global parameters for the simulation. Including
pRun* --Number of orbits to run for (default 500)

pOut* --Frequency with which to save simulation data (default is every 0.1 orbits). 

keep_bins* -- Whether to keep binary stars in the simulation. Deleting any primordial binaries usually speeds up the simulation (default False). 

name -- Prefix for output files (default is archive. Simulation data is ouput to file called archive+tag.bin
is long string of unique letters and numbers).

rt -- Tidal radius for central body (default 1e-4) Useful for detecting collisions. 

coll -- Rebound collision detection algorithm to use (default line)

gravity -- Rebound gravity algorithm to use (default 'basic')

integrator -- Rebound integrator to use (default ias15)

Star indicate the most important parameters...

Following this section you can add an arbitrary number of sections representing different components to add to the 
simulation (e.g. stars, black holes). Note you can call the component anything you want (though it probably would be 
good to pick sensible names). For each section you specify:
N*: # of particles (default 100)

e*: eccentricity (default 0.7)

a_min*: min semi-major axis (sma) (default 1)

a_max*: max sma (default 2)

p: power law slope of sma distribution (default 1--corresponds to an r^-2 surface density profile). 

ang: specifies angle (in degrees) by which j and e vectors (default two degrees). For precisely 
we draw three angles from a normal distribution with standard ang (angle1, angle2, angle3). Then 
rotate jhat (unit vector in direction of angular momentum) by angle1 over major axis and angle minor axis
rotate ehat (unit eccentricity vector) by angle2 over minor axis and angle3 about jhat 

m*: mass of each particle in compoent (default 5e-5)

In addition to these components the code always add a central object with mass 1. 

An example configuration file is provided in the repo (called config). It adds 100 star particles with all of the default
parameters, and a 100 times more massive perturber with 100 times the mass (5e-3). It does not delete the binaries 
and runs for a very short time (10^-8 orbits). 

Output: initial xyz positions, velocities, and masses of all particles in the simulation (in that order) are 
written to init_disk. 

It writes all of the simulation data to a rebound simulation archive (see documentation 
here: https://rebound.readthedocs.io/en/latest/ipython/SimulationArchive.html). Simulation 
archives store all of the simulation data. Afterwards you can load snapshots from an archive and 
calculate all of the orbital elements. 

The code also prints the version of rebound and any detected close encounters with central mass (e.g. TDEs). 







