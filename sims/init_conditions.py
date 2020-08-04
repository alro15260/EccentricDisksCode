import sys
import rebound
import numpy as np
# from astropy.table import Table

name=sys.argv[1]
print name
sa=rebound.SimulationArchive(name)

elem_names=['a', 'e', 'inc', 'omega', 'theta']
#elem_names=['a', 'e']
sim=sa[0]
nn=len(sim.particles)

orbits=sim.calculate_orbits(primary=sim.particles[0])
elems=np.array([getattr(orbits[jj],en) for jj in range(nn-1) for en in elem_names]).reshape([nn-1, len(elem_names)])
# tab=Table(elems, names=elem_names)	
	
# tab.write(name.replace('.bin', '_init.csv'))
np.savetxt(name.replace('.bin', '_init.csv'), elems)