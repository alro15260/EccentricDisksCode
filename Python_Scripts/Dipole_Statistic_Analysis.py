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
import ipyvolume as ipv
import ipyvolume.pylab as p3
import subprocess
from rebound import hash as h
from astropy.io import ascii

%pylab inline
%matplotlib inline

pi=np.pi

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

def get_eccs_dat(sma,ang,simpath):
    '''
    Given simulation dataset, returns the eccentricities of the dataset
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    simpath: base directory location of simulations
    
    Returns:
    Eccentricities of dataset
    '''
    os.chdir(simpath)
    bases=get_bases()
    allmeans=[]
    for j in range(len(bases)):
        datfiles=[]
        for file in glob.glob('*{0}'.format(bases[j])+'_orb_*'):
            datfiles.append(file)
        datfiles.sort(key=natural_keys)
        means=[]
        for i in range(len(datfiles)):
            vals=np.genfromtxt(datfiles[i])[:,1][:-1]
            opp=(vals<1)
            filtered=vals[opp]
            means.append(np.mean(filtered))
        allmeans.append(means)
    mean=np.average(allmeans,axis=0)
    allstds=[]
    allmeans=np.array(allmeans)
    for i in range(2000):
        val=np.std(allmeans[:,i])
        allstds.append(val)
    return mean,allstds

def get_incs_dat(sma,ang,simpath):
    '''
    Given simulation dataset, returns the inclinations of the dataset
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    simpath: base directory location of simulations
    
    Returns:
    Inclinations of dataset
    '''
    os.chdir(simpath)
    bases=get_bases()
    allmeans=[]
    for j in range(len(bases)):
        datfiles=[]
        for file in glob.glob('*{0}'.format(bases[j])+'_orb_*'):
            datfiles.append(file)
        datfiles.sort(key=natural_keys)
        means=[]
        for i in range(len(datfiles)):
            incs=np.genfromtxt(datfiles[i])[:,2][:-1]
            eccs=np.genfromtxt(datfiles[i])[:,1][:-1]
            opp=(eccs<1)
            filtered=incs[opp]
            means.append(np.mean(filtered))
        allmeans.append(means)
    mean=np.average(allmeans,axis=0)
    allstds=[]
    allmeans=np.array(allmeans)
    for i in range(2000):
        val=np.std(allmeans[:,i])
        allstds.append(val)
    return mean*(180/pi),np.array(allstds)*(180/pi)
    def plotter(data,orbit_element,title):
    plt.plot(data[0],color='black')
    x=np.arange(0,2000,1)
    plt.errorbar(x,data[0],data[1],color='red',alpha=0.3)
    plt.xlabel('Orbit Number')
    plt.ylabel('{0}'.format(orbit_element))
    plt.title('{0}'.format(title))
    plt.show()
    
    
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
    Given a 6 column list,returns each of the keplerian orbital elements
    Parameters: A 6 column list referring to each of the 6 keplerian orbital elements
    
    Returns: a tuple(6-D) of each of the orbital elements
    '''
    #given a 6 column list,returns each of the keplerian orbital elements
    a=data[:,0]
    e=data[:,1]
    i=data[:,2]
    Omega=data[:,3]
    omega=data[:,4]
    trueanomaly=data[:,5]
    return a,e,i,Omega,omega,trueanomaly

def get_eccs_finegrid(sma,ang,loc,mass,snap,base):
    '''
    Given simulation parameters, returns the eccentricities of the dataset directly from Rebound
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    mass: Mass of perturber
    snap: Snapshot of simulations out of 500 orbits
    base: base simulation
    
    Returns:
    Eccentricities of dataset
    '''
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
    eccs=[]
    for l,pp in enumerate(p):
        pos=np.array(pp.xyz)-np.array(bh.xyz)
        v=np.array(pp.vxyz)-np.array(bh.vxyz)
        j=np.cross(pos,v)
        rhat=pos/np.linalg.norm(pos)
        e=np.cross(v,j)-rhat
        eccs.append(e)
    return eccs

def enorm(sma,ang,loc,mass,snap,base):
    '''
    Given simulation parameters, returns the normed eccentricity vectors of the dataset directly from Rebound
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    mass: Mass of perturber
    snap: Snapshot of simulations out of 500 orbits
    base: base simulation
    
    Returns:
    Normed Eccentricity vectors of dataset
    '''
    e=get_eccs_finegrid(sma,ang,loc,mass,snap,base)
    allmags=[]
    for i in range(len(e)):
        allmags.append(np.linalg.norm(e[i]))
    allmags=np.array(allmags)
    indices=np.argwhere(allmags>1)
    for i in range(len(indices)):
        #del e[i]
        e=np.delete(e,i,axis=0)
    ehat=[]
    for j in range(len(e)):
        ehat.append(e[j]/(np.linalg.norm(e[j])))
    return np.array(ehat)

def e_par(ehat0):
    '''
    Given unit eccentricity vecotrs, calculates the components in the xy plane
    
    Parameters: Unit eccentricity vectors
    Returns: Unit eccentricity vectors projection in xy plane
    '''
    ehat0[:,2]=0
    return ehat0

def ecc_stat_xy(sma,ang,loc,mass,snap,base):
    '''
    Given simulation parameters, returns the normed eccentricity vectors of the dataset directly from Rebound
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    mass: Mass of perturber
    snap: Snapshot of simulations out of 500 orbits
    base: base simulation
    
    Returns:
    Rayleigh Dipole Statistic in xy plane
    '''
    e_xy=e_par(enorm(sma,ang,loc,mass,snap,base))
    return np.linalg.norm(np.sum(e_xy,axis=0)) / len(e_xy)

def ecc_stat_xy_overbases(sma,ang,loc,mass,snap):
    '''
    Given simulation parameters, returns the normed eccentricity vectors of the dataset directly from Rebound
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    mass: Mass of perturber
    snap: Snapshot of simulations out of 500 orbits
    
    Returns:
    Rayleigh Dipole Statistic in xy plane over all base simulations
    '''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    alleccstat=[]
    for i in range(len(bases)):
        alleccstat.append(ecc_stat_xy(sma,ang,loc,mass,snap,bases[i]))
    return alleccstat

def ecc_stat_xy_all(sma,ang,loc,snap):
    '''
    Given simulation parameters, returns the normed eccentricity vectors of the dataset directly from Rebound
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    snap: Snapshot of simulations out of 500 orbits
    
    Returns:
    Rayleigh Dipole Statistic in xy plane over all base simulations and perturber parameters
    '''
    alleccstats=[]
    alleccstaterrs=[]
    redirs=directories(sma,ang,loc)[0]
    for k in range(len(redirs)):
        mypath=loc+'/'+redirs[k]
        eccstat=ecc_stat_xy_overbases(sma,ang,loc,k,snap)
        alleccstats.append(np.mean(eccstat))
        alleccstaterrs.append(np.std(eccstat))
    return alleccstats,alleccstaterrs

def get_all_data(sma,ang,loc):
    '''
    Given simulation parameters, returns the Rayleigh dipole statistic in the xy plane over all orbits (assumed to be 2000 here)
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    
    Returns:
    Rayleigh Dipole Statistic in xy plane over all orbits
    '''
    data=[]
    for i in range(2000):
        data.append(ecc_stat_xy_all(sma,ang,loc,i))
        print('Orbit {0} complete'.format(i))
    return data

def get_mass_data(data,mass):
    '''
    Given the simulation data, get the Rayleigh Dipole Statistic in the xy plane and uncertainties over masses
    Parameters:
    data: Data from get_all_data function
    mass: Perturber mass to consider
    
    Returns:
    Organized data by mass
    '''
    alldat=[]
    allerr=[]
    for i in range(len(data)):
        dat=data[i][0][mass]
        err=data[i][1][mass]
        alldat.append(dat)
        allerr.append(err)
    return alldat,allerr

def ecc_stat_xy_heavy(sma,ang,loc,snap):
    '''
    Given the simulation data, get the Rayleigh Dipole Statistic in the xy plane and uncertainties for heaviest mass
    Parameters:
    sma: sma or perturber
    ang: inclination of perturber
    loc: base directory location of simulations
    snap: snapshot in simulations
    
    Returns:
    Organized data for heaviest mass perturber
    '''
    alleccstats=[]
    alleccstaterrs=[]
    redirs=directories(sma,ang,loc)[0][-1]
    mypath=loc+'/'+redirs
    eccstat=ecc_stat_xy_overbases(sma,ang,loc,-1,snap)
    alleccstats.append(np.mean(eccstat))
    alleccstaterrs.append(np.std(eccstat))
    return alleccstats,alleccstaterrs

def get_signalstrength(data00):
    f,Pxx=signal.periodogram(data00)
    plt.plot(f,Pxx)
    plt.xlabel('Frequency(1/orbital period)')
    plt.ylabel('Power Spectral Density')
    plt.show()
    print('Period of signal is {:.2f} orbital periods'.format(1/f[np.argmax(Pxx)]))
    
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
    def get_tdename():
    name=glob.glob('*_tde')
    return name

def altcolumnmethod(array,index):
    vals=[]
    if len(array.shape)==2:
        for i in range(len(array)):
            vals.append(array[i][index])
        return vals
    if len(array.shape)==1:
        return array[index] 
        def get_orbital_elements(data):
    #given a 6 column list,returns each of the keplerian orbital elements
    a=data[:,0]
    e=data[:,1]
    i=data[:,2]
    Omega=data[:,3]
    omega=data[:,4]
    trueanomaly=data[:,5]
    return a,e,i,Omega,omega,trueanomaly

def get_orbital_elements_fromrebound(orbits):
    '''
    Given Rebound orbital data,returns each of the keplerian orbital elements in order of a,e,i,Omega,omega, and trueanomaly
    Parameters: Rebound orbital data
    
    Returns: Each of the Keplerian orbital elements'''
    p=[pval.pomega for pval in orbits[:-1]]
    sma=[aval.a for aval in orbits[:-1]]
    eccs=[eccval.e for eccval in orbits[:-1]]
    w=[wval.omega for wval in orbits[:-1]]
    W=[Wval.Omega for Wval in orbits[:-1]]
    return p,sma,eccs,w,W

def get_elements_fromdat(dat):
    '''
    Given a dat file, gets relevant orbital elements
    Paraneters:
    dat: Data file
    
    Returns: relevant orbital elements calculated using Rebound'''
    a,e,i,Omega,omega,trueanomaly=get_orbital_elements(dat)
    sim=rebound.Simulation()
    #Add bh
    sim.add(m = 1)
    #Add 100 stars
    for j in range(len(dat)):
        sim.add(a=a[j],e=e[j],inc=i[j],Omega=Omega[j],omega=omega[j],f=trueanomaly[j])
    orbits=sim.calculate_orbits(primary=sim.particles[0])
    return get_orbital_elements_fromrebound(orbits)

def get_filtered(data,element):
    '''
    Given Rebound orbital data, gets the given orbital element specified by element parameter
    Parameters: 
    data: Rebound orbital data(expecting order of p,sma,e,w,W)
    element: string form of orbital element
    
    Returns: Calculated data for disk of specified element'''
    #expects order of p,sma,e,w,W
    p,sma,e,w,W=get_elements_fromdat(data)
    eccs=np.array(e)
    indices=np.argwhere(eccs>1)
    if element == 'pomega':
        for i in range(len(indices)):
            del p[indices[i][0]]
            p.insert(indices[i][0],np.nan)
        p=np.fmod(np.array(p)+(2*pi),2*pi)
        return p
    if element == 'sma':
        for i in range(len(indices)):
            del sma[indices[i][0]]
            sma.insert(indices[i][0],np.nan)
        return sma
    if element == 'eccentricity':
        for i in range(len(indices)):
            del e[indices[i][0]]
            e.insert(indices[i][0],np.nan)
        return e
    if element == 'pomega2':
        for i in range(len(indices)):
            del w[indices[i][0]]
            w.insert(indices[i][0],np.nan)
        for i in range(len(indices)):
            del W[indices[i][0]]
            W.insert(indices[i][0],np.nan)
        altpomegas=W+w
        altpomegas=np.fmod(np.array(altpomegas)+(2*pi),2*pi)
        return altpomegas

def get_star_orbitalelement(element,sma,ang,loc,mass,snap,base):
    '''
    Gets the orbital element of disk stars under specified perturber parameters
    Parameters:
    element: string form of orbital element
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of simulation data
    mass: mass of perturber
    snap: snapshot of 500 orbits
    base: simulation out of four possible choices
    
    Returns: specified orbital element data for disk under perturber parameters'''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    datfiles=[]
    for file in glob.glob('*{0}'.format(base)+'_orb_*'):
        datfiles.append(file)
    datfiles.sort(key=natural_keys)
    dat=np.genfromtxt(datfiles[snap])
    #order is p,sma,eccs,w,W
    data=get_filtered(dat,element)
    return data
