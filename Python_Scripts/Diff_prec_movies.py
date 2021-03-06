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
import matplotlib.colors as mcol
from scipy import signal
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

#note: various functions include perturber parameters
#due to way data is organized, have random 'perturber' parameters for isolated END so that functions still work...
#... for perturber simulations

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
    #nans not really needed isolated END, since no ejected stars
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


def get_star_orbitalelement_overbases(element,sma,ang,loc,mass,snap):
    '''
    Gets the orbital element of disk stars under specified perturber parameters over all bases
    Parameters:
    element: string form of orbital element
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of simulation data
    mass: mass of perturber
    snap: snapshot of 500 orbits
    
    Returns: specified orbital element data for disk under perturber parameters over all bases'''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    if element == 'pomega':
        dat=[]
        for i in range(len(bases)):
            dat.append(get_star_orbitalelement('pomega',sma,ang,loc,mass,snap,bases[i]))
        return np.nanmean(dat,axis=0)
    if element == 'sma':
        dat=[]
        for i in range(len(bases)):
            dat.append(get_star_orbitalelement('sma',sma,ang,loc,mass,snap,bases[i]))
        return np.nanmean(dat,axis=0)
    if element == 'eccentricity':
        dat=[]
        for i in range(len(bases)):
            dat.append(get_star_orbitalelement('eccentricity',sma,ang,loc,mass,snap,bases[i]))
        return np.nanmean(dat,axis=0)
    

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
    indices=np.argwhere(np.array(allmags)==False)
    indices2=[]
    for i in range(len(indices)):
        indices2.append(indices[i][0])
    indices=indices2
    ie=ie.tolist()
    for i in range(len(indices)):
        del ie[indices[i]]
        ie.insert(indices[i],np.nan)
    return ie
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

#next functions assume data has 5000 orbits
def plotter_snap5000(sma,snap,maxdat0,mindat0,bases0):
    '''
    Plots ie vs sma for a given snapshot while colorcoding the eccentricity of stars
    
    Parameters:
    sma: semi-major axis of perturber
    snap: Snapshot out of the 5000 orbits
    maxdat0: Max eccentricity out of 5000 orbits (calculated elsewhere for colorbar)
    mindat0: Min eccentricity out of 5000 orbits (calculated elsewhere for colorbar)
    bases0: Simulation out of four possible choices
    
    Returns: Plot of ie vs sma at given snapshot'''
    smas=get_star_orbitalelement('sma',sma,180,loc00,-1,snap,bases0[0])
    ies=AngleAnalyzer_data_finegrid(sma,180,loc00,-1,snap,bases0[0])
    eccs=get_star_orbitalelement('eccentricity',sma,180,loc00,-1,snap,bases0[0])
    try:
        c1=plt.scatter(ies,smas,c=eccs,cmap='Blues')
        plt.colorbar(c1,label='e')
    except ValueError:
        pass
    plt.xlabel(r'$i_{e}$')
    plt.ylabel('Semi-Major Axis')
    plt.xlim(-3.4,3.4)
    plt.ylim(0.5,2.5)
    plt.title('Orbit {0}'.format(snap))
    fig=plt.gcf()
    fig.set_size_inches(8.5,6.5)
    
def movie_maker_snap5000(sma,maxdat0,mindat0,bases0):
    '''
    Makes a movie of ie vs sma for a given snapshot while colorcoding the eccentricity of stars
    
    Parameters:
    sma: semi-major axis of perturber
    snap: Snapshot out of the 5000 orbits
    maxdat0: Max eccentricity out of 5000 orbits (calculated elsewhere for colorbar)
    mindat0: Min eccentricity out of 5000 orbits (calculated elsewhere for colorbar)
    bases0: Simulation out of four possible choices
    
    Returns: Movie of ie vs sma'''
    for i in range(5000):
        smas=get_star_orbitalelement('sma',sma,180,loc00,-1,i,bases0[0])
        ies=AngleAnalyzer_data_finegrid(sma,180,loc00,-1,i,bases0[0])
        eccs=get_star_orbitalelement('eccentricity',sma,180,loc00,-1,i,bases0[0])
        try:
            c1=plt.scatter(ies,smas,c=eccs,cmap='Blues')
            plt.colorbar(c1,label='e')
        except ValueError:
            pass
        plt.xlabel(r'$i_{e}$')
        plt.ylabel('Semi-Major Axis')
        plt.xlim(-3.4,3.4)
        plt.ylim(0.5,2.5)
        plt.title('Orbit {0}'.format(i))
        fig=plt.gcf()
        fig.set_size_inches(8.5,6.5)
        plt.savefig('/home/arod/Moviefigs8/Snap{:04d}.png'.format(i))
        plt.clf()
        print('Orbit {0} complete'.format(i))
 
