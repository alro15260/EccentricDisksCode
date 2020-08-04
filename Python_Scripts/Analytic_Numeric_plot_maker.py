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
import uncertainties
from uncertainties import ufloat
%matplotlib inline
pi=np.pi

control0='/home/arod/Research/testclone2/test27_control_manyruns/grid'
control1='/home/arod/Research/testclone2/test27_control_manyruns/grid/M5.0e-02_a2_ang45'

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


#next two functions aren't called

def AngleAnalyzer_finegrid(sma,ang,loc,mass,snap):
    '''
    Function which gets the mean ie values over all four simulations for particular set of pertruber parameters
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    mass: mass of perturber
    snap: Snapshot out of the 500 orbits
    
    Returns: The circular mean of ie for stars in the disk for the four simulations'''
    redirs=directories(sma,ang,loc)[0]
    os.chdir(redirs[mass])
    bases=get_bases()
    circmeans=[]
    for i in range(len(bases)):
        circmeans.append(stats.circmean(AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,bases[i])))
    circmeans=np.array(circmeans)*(180/pi)
    return circmeans

def AngleAnalyzer_all_finegrid(sma,ang,loc,snap):
    '''
    Function which gets the mean ie values for all simulations under different masses
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of perturber simulations
    snap: Snapshot out of the 500 orbits
    
    Returns: The circular mean of ie for stars in the disk for the four simulations for all masses'''
    allmeans=[]
    allerrs=[]
    redirs=directories(sma,ang,loc)[0]
    for k in range(len(redirs)):
        mypath=loc+'/'+redirs[k]
        allcircstds=AngleAnalyzer_finegrid(sma,ang,loc,k,snap)
        allstds.append(np.mean(allcircstds))
        allerrs.append(np.std(allcircstds))
    return allstds,allerrs


def get_orbital_elements(data):
    '''
    Given a 6 column list,returns each of the keplerian orbital elements
    Parameters: A 6 column list referring to each of the 6 keplerian orbital elements
    
    Returns: a tuple(6-D) of each of the orbital elements
    '''
    a=data[:,0]
    e=data[:,1]
    i=data[:,2]
    Omega=data[:,3]
    omega=data[:,4]
    trueanomaly=data[:,5]
    return a,e,i,Omega,omega,trueanomaly

def get_orbital_elements_fromrebound(orbits):
    '''
    Given the orbital data from rebound, gets relevant orbital elements
    Parameters: Orbit data from rebound
    
    Returns:
    p: pomega of stars
    sma: semi-major axis of stars
    eccs: eccentrcity of stars
    w: little omega of stars
    W: big omega of stars'''
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
            #p.insert(indices[i][0],np.nan)
        p=np.fmod(np.array(p)+(2*pi),2*pi)
        return p
    if element == 'sma':
        try:
            for i in range(len(indices)):
                del sma[indices[i][0]]
                #e.insert(indices[i][0],np.nan)
        except:
            IndexError
        return sma
    if element == 'eccentricity':
        try:
            for i in range(len(indices)):
                del e[indices[i][0]]
                #e.insert(indices[i][0],np.nan)
        except:
            IndexError
        return e
    if element == 'pomega2':
        for i in range(len(indices)):
            del w[indices[i][0]]
            #w.insert(indices[i][0],np.nan)
        for i in range(len(indices)):
            del W[indices[i][0]]
            #W.insert(indices[i][0],np.nan)
        altpomegas=W+w
        altpomegas=np.fmod(np.array(altpomegas)+(2*pi),2*pi)
        return altpomegas

#equivilant of get_smas function
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



def get_mean_ies(sma,ang,loc,mass,snap,base,abin):
    '''
    Gets the mean ies of stars in specified semi major axis bins
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of simulation data
    mass: mass of perturber
    snap: snapshot of 500 orbits
    base: simulation out of four possible choices
    abin: semi-major axis range of stars
    
    Returns: averaged ie of stars in specified bins'''
    amin,amax=abin
    ies=AngleAnalyzer_data_finegrid(sma,ang,loc,mass,snap,base)
    smas=get_star_orbitalelement('sma',sma,ang,loc,mass,snap,base)
    smas=np.array(smas)
    return np.mean(ies[(smas>=amin) & (smas<amax)]), np.std(ies[(smas>=amin) & (smas<amax)])

def get_prec_rate(sma,ang,loc,mass,snap,base,abin,step):
    '''
    Gets the numerical precession rate of the disk in a specified bin over some timestep
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of simulation data
    mass: mass of perturber
    snap: snapshot of 500 orbits
    base: simulation out of four possible choices
    abin: semi-major axis range of stars
    step: add step argument which is how far in the future
    
    Returns: disk precession rate'''
    ies1=get_mean_ies(sma,ang,loc,mass,snap,base,abin)
    ies2=get_mean_ies(sma,ang,loc,mass,snap+step,base,abin)
    u=ufloat(ies1).nominal_value
    v=ufloat(ies2).nominal_value
    vec1=[np.cos(u),np.sin(u)]
    vec2=[np.cos(v),np.sin(v)]
    dot=np.dot(vec1,vec2)
    mag=np.linalg.norm(vec1)*np.linalg.norm(vec2)
    theta=np.arccos(dot/mag)
    if theta<2*pi:
        w=(v-u)/step
        return w
    else:
        w=np.fmod((v-u)+(2*pi),2*pi)/step
        return w


def prec_analytic(m_p, a_p, a, e):
    '''
    Analytic prescription (i_e'); note central mass and a of disk particle are assumed to be 1!

    You might want to replace e with a call to a fitting function that you derive by binning numerical data...

    '''
    return (1.5*np.pi*m_p*(a_p/a)**-3.*(1-e**2.)**0.5)*(a**(-3/2))

def get_e_binned(sma,ang,loc,mass,snap,base,abin):
    '''
    Gets the eccentricites of stars within a specified bin
    Parameters:
    sma: semi-major axis of perturber
    ang: inclination of perturber
    loc: base directory of simulation data
    mass: mass of perturber
    snap: snapshot of 500 orbits
    base: simulation out of four possible choices
    abin: semi-major axis range of stars
    
    Returns: the eccentricites of stars within a specified bin
    '''
    amin,amax=abin
    es=get_star_orbitalelement('eccentricity',sma,ang,loc,mass,snap,base)
    smas=get_star_orbitalelement('sma',sma,ang,loc,mass,snap,base)
    es=np.array(es)
    smas=np.array(smas)
    return np.mean(es[(smas>=amin) & (smas<amax)])


def get_analytical_over_sims_over_bins(sma,ang,loc0,mass,snap,bins,a,loc1,ap,mp):
    os.chdir(loc1)
    bases=get_bases()
    all_e_binned=[]
    for i in range(len(bins)):
        dat_over_bases=[]
        for j in range(len(bases)):
            e_binned0=get_e_binned(sma,ang,loc0,mass,snap,bases[j],bins[i])
            dat_over_bases.append(prec_analytic(mp,ap,a[i],e_binned0))
        dat_bin=np.mean(dat_over_bases),(np.std(dat_over_bases)/(np.sqrt(len(bases))))
        all_e_binned.append(dat_bin)
    return all_e_binned

def get_numerical_over_sims_over_bins(sma,ang,loc0,mass,snap,step,bins,loc1):
    os.chdir(loc1)
    bases=get_bases()
    all_prec_rate_binned=[]
    for i in range(len(bins)):
        dat_over_bases=[]
        for j in range(len(bases)):
            prec_rate_binned=get_prec_rate(sma,ang,loc0,mass,snap,bases[j],bins[i],step)
            dat_over_bases.append(prec_rate_binned)
        print('Numerical complete for bins {0}'.format(bins.index(bins[i])))
        dat_bin=np.mean(dat_over_bases),(np.std(dat_over_bases)/(np.sqrt(len(bases))))
        all_prec_rate_binned.append(dat_bin)
    return all_prec_rate_binned

def get_analytical_all(smas,ang,loc0,mass,snap,bins,a,ap,num_data):
    alldat=[]
    m=np.geomspace(5e-2,5e-1,10)
    for i in range(len(smas)):
        redirs=directories(smas[i],ang,loc0)[0]
        dat=(get_analytical_over_sims_over_bins(smas[i],ang,loc0,mass,snap,bins,a,loc00+'/'+redirs[mass],smas[i],m[mass]))
        alldat.append(dat)
    return alldat

def get_num_analytical_flatness(analytical_data,num_data):
    dat=[]
    analytical_data=np.array(analytical_data)
    for i in range(len(analytical_data)):
        dat.append(analytical_data[i][:,0])
    dat=np.array(dat).tolist()
    newdat=[]
    for i in range(len(smas)):
        newdat.append((dat[i]+num_data)[-1]/((dat[i]+num_data)[0]))
    return np.absolute(np.array(newdat)-1)
    #return (np.array(newdat)-1)

def get_num_analytical_flatness_over_masses(smas,ang,loc0,snap,bins,a,ap,num_data):
    allflattests=[]
    m=np.geomspace(5e-2,5e-1,10)
    for i in range(10):
        analytical_dat=get_analytical_all(smas,ang,loc0,i,snap,bins,a,ap,num_data)
        flattest=get_num_analytical_flatness(analytical_dat,num_data)
        allflattests.append(flattest)
    return allflattests

def GridPlotter_num_analytical(smas,ang,loc0,snap,bins,a,ap,num_data):
    masses=[0,1,2,3,4,5,6,7,8,9]
    x=masses
    y=smas
    delta_x=np.diff(x)[0]
    delta_y=np.diff(y)[0]
    dat=get_num_analytical_flatness_over_masses(smas,ang,loc0,snap,bins,a,ap,num_data)
    dat=np.array(dat).tolist()
    dat=np.transpose(dat)
    #for i in range(len(smas)):
        #dat.append(np.absolute(get_num_analytical_flatness_over_masses(smas[i],ang,loc0,snap,bins,a,ap,num_data)))
        #dat.append(get_slope_overmasses(smas[i],ang,loc0,snap,bins,a,smas[i],num_data))
    #dat=np.array(dat)
    fig,ax=plt.subplots(figsize=(10,9))
    #for ii,yy in enumerate(y):
        #for jj,xx in enumerate(x):
            #plt.text(xx, yy, '{0:.0f}'.format(dat[ii, jj]), fontsize=10, horizontalalignment='center')
    x_ords=np.append(x-delta_x/2.0,x[-1]+delta_x/2.0)
    y_ords=np.append(y-delta_y/2.0,y[-1]+delta_y/2.0)
    cmap=cm.get_cmap('Blues')
    cmap.set_under('white')
 
    fulldat=dat
    min0=np.min(fulldat)
    max0=np.max(fulldat)
    ax0=plt.gca()
    p1=ax.pcolormesh(x_ords, y_ords,(dat), cmap=cmap, norm=colors.Normalize(vmin=min0, vmax=max0))
    cbar=fig.colorbar(p1, ax=ax, label=r"$|(i_{e}'(2)/i_{e}'(1))-1|$", extend='both')
    plt.xticks(())
    plt.text(-.75,5.1,'$5 x 10^{-2}$',fontsize=22)
    plt.text(4.2,5.1,'$5 x 10^{-1.5}$',fontsize=22)
    plt.text(9,5.1,'$5 x 10^{-1}$',fontsize=22)
    label=ax0.set_xlabel(r'Perturber Mass $[M_p]$',fontsize=22)
    ax0.xaxis.set_label_coords(.5,-.075)
    plt.ylabel('Semi-Major Axis')
    plt.savefig('/home/arod/PaperPlots/AnalyticalNumerical.pdf')
    plt.show()
