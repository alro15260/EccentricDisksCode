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
pi=np.pi
loc1='/home/arod/Research/testclone2/test19_highmasshighsma/grid2/grid/M5.0e-01_a9.0_ang180'
loc3='/home/arod/Research/testclone2/test19_highmasshighsma/grid2/grid/'

def directories(a,ang,base):
    'Given an angle, returns a sorted list of directories from lowest to highest Mass'
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
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def get_bases():
    fs=glob.glob("*orb_0.dat")
    fs=[ff.replace("_orb_0.dat", "") for ff in fs]
    return fs

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
    p=[pval.pomega for pval in orbits[:-1]]
    sma=[aval.a for aval in orbits[:-1]]
    eccs=[eccval.e for eccval in orbits[:-1]]
    w=[wval.omega for wval in orbits[:-1]]
    W=[Wval.Omega for Wval in orbits[:-1]]
    return p,sma,eccs,w,W

def get_elements_fromdat(dat):
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
    
def get_colors(length):
    x=np.linspace(0,1,length)[::-1]
    y=np.linspace(0,0,length)
    z=np.linspace(0,1,length)
    colors=[]
    for i in range(len(x)):
        color=x[i],y[i],z[i]
        colors.append(color)
    return colors

def AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,base):
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
    
def plotter_snap5000(sma,snap,maxdat0,mindat0,bases0):
    smas=get_star_orbitalelement('sma',sma,180,loc00,-1,snap,bases0[0])
    ies=AngleAnalyzer_data_finegrid(sma,180,loc00,-1,snap,bases0[0])
    eccs=get_star_orbitalelement('eccentricity',sma,180,loc00,-1,snap,bases0[0])
    z=np.linspace(0,1,100)
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
    for i in range(5000):
        smas=get_star_orbitalelement('sma',sma,180,loc00,-1,i,bases0[0])
        ies=AngleAnalyzer_data_finegrid(sma,180,loc00,-1,i,bases0[0])
        eccs=get_star_orbitalelement('eccentricity',sma,180,loc00,-1,i,bases0[0])
        z=np.linspace(0,1,100)
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
    
def plotter_snap4500(sma,snap,maxdat0,mindat0,bases0):
    smas=get_star_orbitalelement('sma',sma,45,loc03,-1,snap,bases0[0])
    ies=AngleAnalyzer_data_finegrid(sma,45,loc03,-1,snap,bases0[0])
    eccs=get_star_orbitalelement('eccentricity',sma,45,loc03,-1,snap,bases0[0])
    z=np.linspace(0,1,99)
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
    
def movie_maker_snap4500(sma,maxdat0,mindat0,bases0):
    for i in range(4500):
        smas=get_star_orbitalelement('sma',sma,45,loc03,-1,i,bases0[0])
        ies=AngleAnalyzer_data_finegrid(sma,45,loc03,-1,i,bases0[0])
        eccs=get_star_orbitalelement('eccentricity',sma,45,loc03,-1,i,bases0[0])
        z=np.linspace(0,1,99)
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
        plt.savefig('/home/arod/Moviefigs9/Snap{:04d}.png'.format(i))
        plt.clf()
        print('Orbit {0} complete'.format(i))
