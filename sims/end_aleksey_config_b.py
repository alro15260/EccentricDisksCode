#!/usr/bin/env python

import ConfigParser
import argparse
import uuid
# import sys

import numpy as np
from collections import OrderedDict
# sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import random as rand
from bin_analysis import bin_find_sim
import math

import os

def rotate_vec(angle,axis,vec):    
	'''
	Rotate vector vec by angle around axis (couter-clockwise)
	'''
	vRot = vec*math.cos(angle) + np.cross(axis,vec)*math.sin(angle) + axis*np.dot(axis,vec)*(1 -math.cos(angle))
	return vRot	

def gen_disk(ang1, ang1_mean, ang2, ang2_mean, ang3, ang3_mean):
	'''
	This is from some old code that starts with perfectly aligned e and j vectors and then rotates them by a small amount
	'''
	ehat = np.array([1,0,0])
	jhat = np.array([0,0,1])
	bhat = np.cross(jhat,ehat)    # rotate jhat by angle1 over major axis and angle 2 minor axis
	# rotate ehat by angle2 over minor axis (for consistency) and angle3 about jhat
	angle1 = np.random.normal(ang1_mean, ang1, 1)
	angle2 = np.random.normal(ang2_mean, ang2, 1)
	angle3 = np.random.normal(ang3_mean, ang3, 1)    
	jhat = rotate_vec(angle1,ehat,jhat)
	jhat = rotate_vec(angle2,bhat,jhat)
	ehat = rotate_vec(angle2,bhat,ehat)
	ehat = rotate_vec(angle3,jhat,ehat)    
	n = np.cross(np.array([0,0,1]), jhat)
	n = n / np.linalg.norm(n)   
	Omega = math.atan2(n[1], n[0])
	omega = math.acos(np.dot(n, ehat))
	if ehat[2] < 0:
		omega = 2*np.pi - omega    
	inc=math.acos(jhat[2])    
	return inc, Omega, omega


def density(min1, max1, p):
	'''
	Generate a random from a truncated power law PDF with power law index p. 
	min1 and max1

	'''
	r=np.random.random(1)[0]
	if p==1:
		return min1*np.exp(r*np.log(max1/min1))
	else:
		return (r*(max1**(1.-p)-min1**(1.-p))+min1**(1.-p))**(1./(1-p))

def heartbeat(sim):
	print(sim.contents.dt, sim.contents.t)
# sim is a pointer to the simulation object,
# thus use contents to access object data.
# See ctypes documentation for details.
	# print(sim.contents.dt)

def get_tde(sim, reb_coll):
	orbits = sim[0].calculate_orbits(primary=sim[0].particles[0])
	p1,p2 = reb_coll.p1, reb_coll.p2
	idx, idx0 = max(p1, p2), min(p1, p2)
	if idx0==0:
		##idx decremented by 1 because there is no orbit 0
		name=sim[0].simulationarchive_filename
		f=open(name.replace('.bin', '_tde'), 'a+')
		f.write('{0} {1} {2} {3} TDE!\n'.format(sim[0].t, orbits[idx-1].a, orbits[idx-1].e, idx))
		f.close()

	return 0

def main():
	parser=argparse.ArgumentParser(
		description='Set up a rebound run')
	parser.add_argument('--config', nargs=1, default='config',
		help='File containing simulation parameters')
	# parser.add_argument('--keep_bins', action='store_true',
	# 	help="Don't delete bins from simulation")


	##Parsing command line arguments.
	args=parser.parse_args()
	config_file=args.config
	##Unique tag for output file.
	tag=str(uuid.uuid4())

	##Default stellar parameters 
	config=ConfigParser.SafeConfigParser(defaults={'name': 'archive'.format(tag), 'N':'100', 'e':'0.7',
		'gravity':'basic', 'integrator':'ias15', 'dt':'0', \
		'a_min':'1.', 'a_max':'2.', 'ang1_mean':'0', 'ang2_mean':'0', 'ang3_mean':'0', 'ang1':'2.', 'ang2':'2.', 'ang3':'2.', 'm':'5e-5', 'keep_bins':'False', 'rt':'1.0e-4', 'coll':'line', 'pRun':'500', 'pOut':'50', 
		'p':'1', 'frac':'nan', 'outDir':'./'}, dict_type=OrderedDict)
	# config.optionxform=str
	config.read(config_file)

	##Name of our put file 
	name=config.get('params', 'name')
	name=name+"_"+tag+".bin"
	name=config.get('params', 'outDir')+'/'+name
	##Length of simulation and interval between snapshots
	pRun=config.getfloat('params', 'pRun')
	pOut=config.getfloat('params', 'pOut')
	keep_bins=config.getboolean('params', 'keep_bins')
	rt=config.getfloat('params', 'rt')
	coll=config.get('params', 'coll')

	print pRun, pOut, rt, coll
	sections=config.sections()
	sections=sections[1:]
	##Initialized the rebound simulation
	sim = rebound.Simulation()
	sim.G = 1.	
	##Central object
	sim.add(m = 1, r=rt) 
	sim.gravity=config.get('params', 'gravity')
	sim.integrator=config.get('params', 'integrator')
	dt=config.getfloat('params', 'dt')
	if dt:
		sim.dt=dt
	if sim.gravity=='tree':
		##Fixing box, angle, and boundary parameters in the tree code.
		sim.configure_box(10.)
		sim.boundary='open'
		sim.opening_angle2=1.5

	buff=1.5
	nparts={}
	num={}
	##Add particles; Can have different sections with different types of particles (e.g. heavy and light)
	##see the example config file in repository. Only require section is params which defines global parameters 
	##for the simulation (pRun and pOut).
	for ss in sections:
		num[ss]=int(config.get(ss, 'N'))
		N=int(buff*num[ss])
		e=config.getfloat(ss, 'e')
		m=config.getfloat(ss, 'm')
		frac=config.getfloat(ss, 'frac')
		if not np.isnan(frac):
			m=frac/num[ss]
		a_min=config.getfloat(ss, 'a_min')
		a_max=config.getfloat(ss, 'a_max')
		p=config.getfloat(ss, 'p')
		ang1_mean=config.getfloat(ss, 'ang1_mean')
		ang1=config.getfloat(ss, 'ang1')
		ang2_mean=config.getfloat(ss, 'ang2_mean')
		ang2=config.getfloat(ss, 'ang2')
		ang3_mean=config.getfloat(ss, 'ang3_mean')
		ang3=config.getfloat(ss, 'ang3')

		N0=len(sim.particles)
		for l in range(0,N): # Adds stars
			##Use AM's code to generate disk with aligned eccentricity vectors, but a small scatter in i and both omegas...
			inc, Omega, omega=gen_disk(ang1*np.pi/180., ang1_mean*np.pi/180., ang2*np.pi/180., ang2_mean*np.pi/180., ang3*np.pi/180., ang3_mean*np.pi/180.0)
			a0=density(a_min, a_max, p)
			M = rand.uniform(0., 2.*np.pi)
			sim.add(m = m, a = a0, e = e, inc=inc, Omega = Omega, omega = omega, M = M, primary=sim.particles[0])
		##Indices of each component
		nparts[ss]=(N0,N0+N-1)
	
	f=open('init_disk', 'w')
	sim.move_to_com()
	for ii in range(len(sim.particles)):
		f.write('{0:.16e} {1:.16e} {2:.16e} {3:.16e} {4:.16e} {5:.16e} {6:.16e}\n'.format(sim.particles[ii].x, sim.particles[ii].y, sim.particles[ii].z,\
			sim.particles[ii].vx, sim.particles[ii].vy, sim.particles[ii].vz, sim.particles[ii].m))
	f.close()

	fen=open(name.replace('.bin', '_en'), 'a')
	fen.write(sim.gravity+'_'+sim.integrator+'_'+'{0}'.format(sim.dt))
	if not keep_bins:
		##Integrate forward a small amount time to initialize accelerations.
		sim.move_to_com()
		sim.integrate(1.0e-15)
		##Look for binaries
		bins=bin_find_sim(sim)
		bins=np.array(bins)
		#print len(bins[:,[1,2]])
		##Delete all the binaries that we found. The identification of binaries depends in part on the tidal field 
		##of the star cluster, and this will change as we delete stars. So we repeat the binary 
		##deletion process several times until there are none left.
		while len(bins>0):
			##Delete in reverse order (else the indices would become messed up)
			to_del=(np.sort(np.unique(bins[:,1]))[::-1]).astype(int)
			print "deleting",len(to_del)
			for idx in to_del:
				sim.remove(idx)
			sim.integrate(sim.t+sim.t*1.0e-14)
			bins=bin_find_sim(sim)
			N0=1
			##Update indices for each section after binary deletion
			for ss in sections:
				del1=len(np.intersect1d(range(nparts[ss][0],nparts[ss][-1]+1), to_del))
				tot1=nparts[ss][-1]-nparts[ss][0]+1
				nparts[ss]=(N0, N0+tot1-del1-1)
				N0=N0+tot1-del1

	##Delete all of the excess particles
	for ss in sections[::-1]:
		to_del=range(nparts[ss][0]+num[ss], nparts[ss][-1]+1)[::-1]
		for idx in to_del:
			sim.remove(idx)
	print len(sim.particles)

	ms=np.array([pp.m for pp in sim.particles])
	print ms
	sim.collision=coll
	sim.collision_resolve=get_tde

	##Set up simulation archive for output
	sim.automateSimulationArchive(name,interval=2.0*np.pi*pOut,deletefile=True)
	#sim.heartbeat=heartbeat
	sim.move_to_com()
	sim.simulationarchive_snapshot(name)


	en=sim.calculate_energy()
	print sim.N, rebound.__version__
	t=0.0
	orb_idx=0
	while(t<pRun*2.0*np.pi):
		orbits=sim.calculate_orbits(primary=sim.particles[0])
		np.savetxt(name.replace('.bin', '_orb_{0}.dat'.format(orb_idx)), [[oo.a, oo.e, oo.inc, oo.Omega, oo.omega, oo.f] for oo in orbits])
		sim.integrate(sim.t+2.0*np.pi)
		t+=2.0*np.pi
		orb_idx+=1
	# sim.integrate(pRun*2*np.pi)
	en2=sim.calculate_energy()
	print abs(en2-en)/en
	fen.write('_{0:2.3g}'.format(abs(en2-en)/en))
	fen.close()





if __name__ == '__main__':
	main()








