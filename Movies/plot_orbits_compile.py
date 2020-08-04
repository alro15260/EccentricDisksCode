import rebound 
import numpy as np
import sys
import shlex

from PyAstronomy import pyasl
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import subprocess
import rebound

def bash_command(cmd):
	'''Run command from the bash shell'''
	process=subprocess.Popen(['/bin/bash', '-c',cmd],  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	return process.communicate()[0]

##Location of simulation data...
loc=sys.argv[1]
##You will have to modify the names of the data files (the archive_out part) here. 
frames=bash_command('echo {0}/archive_out_*.dat'.format(loc)).decode('utf-8')
frames=len(shlex.split(frames))
##and here
samp=np.genfromtxt(loc+'/archive_out_0.dat')

pts=200
x=np.empty([frames, len(samp), pts+1])
y=np.empty([frames, len(samp), pts+1])
z=np.empty([frames, len(samp), pts+1])
sim=rebound.Simulation()
sim.add(m=4e6)

for idx in range(frames):
	##and here.
	dat=np.genfromtxt(loc+'/archive_out_{0}.dat'.format(idx))

	for i in range(len(dat)):
		if dat[i,0]<0:
			continue
		sim.add(a=dat[i,0], e=dat[i,1], inc=dat[i,2], Omega=dat[i, 3], omega=dat[i,4], f=dat[i, 5], m=100)
		p=sim.particles[1]
		pos=np.array(p.sample_orbit(pts, primary=sim.particles[0]))
		sim.remove(1)

		x[idx,i,:-1] = pos[:,0]
		x[idx,i][-1]=np.nan
		y[idx,i,:-1] = pos[:,1]
		y[idx,i][-1]=np.nan
		z[idx,i,:-1] = pos[:,2]
		z[idx,i][-1]=np.nan
np.savez(loc+'/x.npz', x)
np.savez(loc+'/y.npz', y)
np.savez(loc+'/z.npz', z)
