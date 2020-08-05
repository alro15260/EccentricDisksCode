import numpy as np
import os
import scipy
import matplotlib.pyplot as plt
from scipy import stats
import rebound
import glob
from matplotlib import ticker,cm
from matplotlib import axis
import matplotlib.colors as colors
from collections import OrderedDict
import re
from scipy import signal
%matplotlib inline

%pylab inline
%matplotlib inline
pi=np.pi

#Defining Necessary Functions

#Sort Directories by Mass function
def directories(a,ang,base):
    '''
    Given parameters of a given simulation, returns a list of directories sorted by mass(lowest to highest)
    Parameters:
    a: sma or perturber
    ang: inclination of perturber
    base: base directory location of simulations
    
    Returns:
    List of directories sorted from lowest to highest perturber mass
    '''
    dirs=[]
    os.chdir(base)
    for file in glob.glob('*_a{0}_ang{1}'.format(a,ang)):
        dirs.append(file)
    vals=[]
    for i in range(len(dirs)):
        val=float(dirs[i][1:8])
        vals.append(val)
    newvals=np.sort(np.array(vals))
    newvals=list(newvals)
    newvals0=[]
    for i in range(len(newvals)):
        newvals0.append('{:.1e}'.format(newvals[i]))
    redirs=[]
    for i in range(len(dirs)):
        redirs.append('M'+newvals0[i]+dirs[i][8:])
    return redirs,newvals

def get_bases():
    '''
    Gets the names of the four different simulation names of one particular set of parameters
    Parameters: None
    
    Returns: List of four different simulation names '''
    fs=glob.glob("*orb_0.dat")
    fs=[ff.replace("_orb_0.dat", "") for ff in fs]
    return fs

def atoi(text):
    '''
    Function which turns text to integers to be more easily sorted
    Parameters: Text to be sorted
    
    Returns: integer datatype of text'''
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    Sorts text from lowest to highest so dat files are more easily read
    Parameters: data file in text format
    
    Returns: List of sorted datfiles
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def get_bases():
    '''
    Gets the names of the four different simulation names of one particular set of parameters
    Parameters: None
    
    Returns: List of four different simulation names '''
    fs=glob.glob("*orb_0.dat")
    fs=[ff.replace("_orb_0.dat", "") for ff in fs]
    return fs

def get_orbital_elements(data):
    '''
    Given a 6 column list,returns each of the keplerian orbital elements in order of a,e,i,Omega,omega, and trueanomaly
    Parameters: 6-D list to analyze
    
    Returns: Each of the Keplerian orbital elements'''
    a=data[:,0]
    e=data[:,1]
    i=data[:,2]
    Omega=data[:,3]
    omega=data[:,4]
    trueanomaly=data[:,5]
    return a,e,i,Omega,omega,trueanomaly

def AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,base):
    '''
    For particular perturber parameters, simulation, and snapshot, returns ie values of the disk
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    mass: mass of perturber
    snap: Snapshot out of the 500 orbits
    base: Simulation out of four possible choices
    
    Returns: ie values for the disk'''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    datfiles=[]
    for file in glob.glob('*{0}'.format(base)+'_orb_*'):
        datfiles.append(file)
    datfiles.sort(key=natural_keys)
    dat=np.genfromtxt(datfiles[snap])
    a,e,i,Omega,omega,trueanomaly=get_orbital_elements(dat)
    sim=rebound.Simulation()
    #Add bh
    sim.add(m = 1)
    #Add 100 stars
    for j in range(len(dat)):
        sim.add(a=a[j],e=e[j],inc=i[j],Omega=Omega[j],omega=omega[j],f=trueanomaly[j])
    p=sim.particles[1:-1]
    bh=sim.particles[0]
    ie=np.zeros(len(p))
    allmags=[]
    #initialize array
    for l,pp in enumerate(p):
        pos=np.array(pp.xyz)-np.array(bh.xyz)
        v=np.array(pp.vxyz)-np.array(bh.vxyz)
        j=np.cross(pos,v)
        rhat=pos/np.linalg.norm(pos)
        e=np.cross(v,j)-rhat
        mag=np.linalg.norm(e)
        allmags.append(mag<1)
        ie[l]=np.arctan2(e[1],e[0])
    return ie[allmags]

def AngleAnalyzer_finegrid(sma,ang,loc,mass,snap):
    '''
    For particular perturber parameters, simulation, and snapshot, returns ie values of the disk over base simulations
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    mass: mass of perturber
    snap: Snapshot out of the 500 orbits
    
    Returns: ie values for the disk over base simulations'''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    circstds=[]
    for i in range(len(bases)):
        circstds.append(stats.circstd(AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,bases[i])))
    circstds=np.array(circstds)*(180/pi)
    return circstds

def AngleAnalyzer_all_finegrid(sma,ang,loc,snap):
    '''
    For particular perturber parameters, simulation, and snapshot, returns ie values of the disk over base simulations and masses
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    snap: Snapshot out of the 500 orbits
    
    Returns: ie values for the disk over base simulations over perturber masses'''
    allstds=[]
    allerrs=[]
    redirs=directories(sma,ang,loc)[0]
    for k in range(len(redirs)):
        mypath=loc+'/'+redirs[k]
        allcircstds=AngleAnalyzer_finegrid(sma,ang,loc,k,snap)
        allstds.append(np.mean(allcircstds))
        allerrs.append(np.std(allcircstds))
    return allstds,allerrs

def get_dat(sma,ang,loc):
    '''
    For particular perturber parameters, simulation, and snapshot, returns ie values of the disk over snapshots (assumed to be 2000 here)
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    
    Returns: ie values for the disk over base simulations over perturber masses over all orbits'''
    allvals=[]
    allerrs=[]
    for i in range(2000):
        allvals.append(AngleAnalyzer_all_finegrid(sma,ang,loc,i)[0][0])
        allerrs.append(AngleAnalyzer_all_finegrid(sma,ang,loc,i)[1][0])
        print('Orbit {0} complete'.format(i))
    return allvals,allerrs
