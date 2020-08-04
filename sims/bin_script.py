from rebound_runs import bin_analysis
import rebound
import sys
import numpy as np


from astropy.table import Table

def vh(m, a):
	return m**(1./3.)*a**-0.5

name=sys.argv[1]
# bins=bin_analysis.BinAnalysis(name)
##Set up simulation archive
#--------------------------------------------------------------------------------------------------_#
sa=rebound.SimulationArchive(name)
ms = np.array([pp.m for pp in sa[0].particles[1:]])
#ms=np.genfromtxt(name.replace('.bin', '_masses'))
#--------------------------------------------------------------------------------------------------_#
##Velocity dispersions to file
sigs=np.empty([len(sa), 3])
sigs_low=np.empty([len(sa), 3])
sigs_high=np.empty([len(sa), 3])

##Save orbital elements to hdf5 file.
#--------------------------------------------------------------------------------------------------_#
elem_names=['a', 'e', 'inc', 'omega']
#elem_names=['a', 'e']
interval=1
# sigs_low_b=np.empty(len(sa))

for ii in range(0, len(sa), interval):
	sim=sa[ii]
	vs = np.array([np.array(pp.vxyz) for pp in sim.particles[1:]])

	nn=len(sim.particles)

	orbits=sim.calculate_orbits(primary=sim.particles[0])
	elems=np.array([getattr(orbits[jj],en) for jj in range(nn-1) for en in elem_names]).reshape([nn-1, len(elem_names)])
	tab=Table(elems, names=elem_names)	
	
	tab.write(name.replace('.bin','_elems.hdf5'), '/{0}'.format(ii), format='hdf5', append=True, overwrite=True)
	##Low mass stars
	filt1=(ms<=np.median(ms))
	filt1b=(ms>np.median(ms))
	##Ignore any unbound stars
	filt2=(tab['a']>0)
	#sigs_low_b[ii]=np.std(vs[filt1 & filt2][:,2]/vh(np.min(ms), tab['a'][filt1 & filt2]))
	sigs[ii] = np.std(vs[filt2],axis=0) 
	sigs_high[ii] = np.std(vs[filt1b & filt2], axis=0)
	sigs_low[ii] = np.std(vs[filt1 & filt2], axis=0)


np.savetxt(name.replace('.bin', '_sigs'), sigs)
np.savetxt(name.replace('.bin', '_sigs_high'), sigs_high)
np.savetxt(name.replace('.bin', '_sigs_low'), sigs_low)
#np.savetxt(name.replace('.bin', '_vvh_ratio_low'), sigs_low_b)