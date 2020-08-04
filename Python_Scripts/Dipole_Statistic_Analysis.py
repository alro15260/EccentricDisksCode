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

def get_eccs_dat(sma,ang,simpath):
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

def get_eccs_finegrid(sma,ang,loc,mass,snap,base):
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
    ehat0[:,2]=0
    return ehat0

def ecc_stat_xy(sma,ang,loc,mass,snap,base):
    e_xy=e_par(enorm(sma,ang,loc,mass,snap,base))
    return np.linalg.norm(np.sum(e_xy,axis=0)) / len(e_xy)

def ecc_stat_xy_overbases(sma,ang,loc,mass,snap):
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    alleccstat=[]
    for i in range(len(bases)):
        alleccstat.append(ecc_stat_xy(sma,ang,loc,mass,snap,bases[i]))
    return alleccstat

def ecc_stat_xy_all(sma,ang,loc,snap):
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
    data=[]
    for i in range(2000):
        data.append(ecc_stat_xy_all(sma,ang,loc,i))
        print('Orbit {0} complete'.format(i))
    return data

def get_mass_data(data,mass):
    alldat=[]
    allerr=[]
    for i in range(len(data)):
        dat=data[i][0][mass]
        err=data[i][1][mass]
        alldat.append(dat)
        allerr.append(err)
    return alldat,allerr

def ecc_stat_xy_heavy(sma,ang,loc,snap):
    alleccstats=[]
    alleccstaterrs=[]
    redirs=directories(sma,ang,loc)[0][-1]
    mypath=loc+'/'+redirs
    eccstat=ecc_stat_xy_overbases(sma,ang,loc,-1,snap)
    alleccstats.append(np.mean(eccstat))
    alleccstaterrs.append(np.std(eccstat))
    return alleccstats,alleccstaterrs

def get_heavy_data(sma,ang,loc):
    data=[]
    keys=[1,2,3,4,5]
    for i in range(2000):
        data.append(ecc_stat_xy_heavy(sma,ang,loc,i))
        print('Orbit {0} complete'.format(i))
    return data
def t_osc(e,sigma_ie):
    t_sec=200
    return ((e*sigma_ie)**(3/2))*(t_sec)


def get_eccs_dat2(sma,ang,simpath,mass):
    redirs=directories(sma,ang,simpath)[0]
    os.chdir(redirs[mass])
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

def get_all_tosc(sma,ang,loc,mass):
    e=get_eccs_dat2(sma,ang,loc,mass)
    t_osc_all=[]
    t_osc_err=[]
    for i in range(2000):
        eorbit=e[0][i]
        eorbit_err=e[1][i]
        sigma_ie,sigma_ie_err=np.array(AngleAnalyzer_all_finegrid(sma,ang,loc,i))[:,mass]*(pi/180)
        u=ufloat(sigma_ie,sigma_ie_err)
        v=ufloat(eorbit,eorbit_err)
        w=t_osc(u,v)
        wval=w.nominal_value
        werr=w.std_dev
        t_osc_all.append(wval)
        t_osc_err.append(werr)
        print('Orbit {0} complete'.format(i))
    return t_osc_all,t_osc_err

def get_signalstrength(data00):
    f,Pxx=signal.periodogram(data00)
    plt.plot(f,Pxx)
    plt.xlabel('Frequency(1/orbital period)')
    plt.ylabel('Power Spectral Density')
    plt.show()
    print('Period of signal is {:.2f} orbital periods'.format(1/f[np.argmax(Pxx)]))
    
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
    return ie[allmags]

def AngleAnalyzer_finegrid(sma,ang,loc,mass,snap):
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    circstds=[]
    for i in range(len(bases)):
        circstds.append(stats.circstd(AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,bases[i])))
    circstds=np.array(circstds)*(180/pi)
    return circstds

def AngleAnalyzer_all_finegrid(sma,ang,loc,snap):
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

def get_inneredge(sma,ang,loc,mass,snap,base):
    smas=get_star_orbitalelement('sma',2,45,loc1,mass,snap,bases[0])
    pomegas=get_star_orbitalelement('pomega',2,45,loc1,mass,snap,bases[0])
    indices=np.argsort(smas)[:76]
    newpomegas=[]
    for i in range(len(indices)):
        newpomegas.append(pomegas[indices[i]])
    return newpomegas

def completion_per(value):
    pers=np.arange(0,2100,200)
    if value in pers:
        print('{:.0f}% complete'.format((value/2000)*100))
        
def get_all_orbits(element,index):
    os.chdir(loc2)
    bases=get_bases()
    allorbits=[]
    for i in range(2000):
        allorbits.append(get_star_orbitalelement(element,2,45,loc1,-1,i,bases[0])[int(index)])
        completion_per(i)
    return allorbits

def get_all_orbits_inneredge():
    os.chdir(loc2)
    bases=get_bases()
    allorbits=[]
    for i in range(2000):
        allorbits.append(np.mean(get_inneredge(2,45,loc1,-1,i,bases[0])))
        completion_per(i)
    return allorbits

def get_pomega_data(loc,index,edge):
    os.chdir(loc)
    names=get_tdename()
    bases=get_bases()
    bases.sort(key=natural_keys)
    names.sort(key=natural_keys)
    indices=[]
    #use only first tde file for now
    dat=np.genfromtxt(names[0])
    indices.append(np.unique(dat[:,-2]))
    dat=[]
    for j in range(len(indices)):
        dat.append(indices[j].tolist())
    dat=dat[0]
    orbits=get_all_orbits('pomega',index)
    return np.array(orbits)-np.array(edge)
    
 def tosc_plotter(indices):
    inneredge=get_all_orbits_inneredge()
    for i in range(len(indices)):
        e=get_all_orbits('eccentricity',indices[i])
        p=get_pomega_data(loc2,indices[i],inneredge)
        p0=p-inneredge
        dat=[]
        for j in range(2000):
            orbit=p0[j],e[j]
            dat.append(tuple(orbit))
        dat=np.array(dat)
        plt.scatter(dat[:,0]*(180/pi),dat[:,1])
        plt.xlabel(r'$\bar{\omega}-<\bar{\omega}_{inner\,\,disk}>(deg)$')
        plt.ylabel('e')
        plt.show()
