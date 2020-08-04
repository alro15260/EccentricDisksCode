import numpy as np
import os
import scipy
import matplotlib.pyplot as plt
from scipy import stats
import rebound
import glob
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

mypath2='/home/arod/Research/testclone2/test13_control/grid/grid_run12/grid_run12/grid/'
# os.chdir(mypath)
#Get random generated string coorseponding to each trial
def get_bases():
    fs=glob.glob("*orb_0.dat")
    fs=[ff.replace("_orb_0.dat", "") for ff in fs]
    return fs

#Given the name of a subdirectory within the mypath directory, unpack the data into an array of shape: (4,101,500,6)
def get_stacked_data(name):
    cwd = os.chdir(name)
    bases = get_bases()
    print(bases)

    dat_stacked=np.empty([len(bases), 100, 500, 6])

    for i1,base in enumerate(bases):
        for i2 in range(500):
            try:
                print(base+'_orb_{0}'.format(i2)+'.dat')
                dat_stacked[i1, :,i2]=np.genfromtxt(base+'_orb_{0}'.format(i2)+'.dat')
            except Exception as e:
                print(e)
                continue
    return dat_stacked

def get_eccs3(sma,ang,base):
    ##Getting data for all trials and snapshots at a particular location.
    #Highest mass
    loc=directories(sma,ang,base)[0][0]

    print(base+'/'+loc)
    dat=get_stacked_data(base+'/'+loc)
    ##Eccentricity--1st index is trial number, second index is particle number, third index is snapshot
    ##and final index is the orbital element.
    eccs=dat[:,:,:,1]
    ##Number of trials and snapshots
    trials=eccs.shape[0]
    snapshots=eccs.shape[-1]
    ##Store the mean for eccentricity for all trials and snapshots
    mean=np.ones([trials, snapshots])*np.nan
    for i in range(trials):
        ##Get ith trial
        tmp=eccs[i]
        print(eccs)
        ##Flip indices--snapshot number then particle
        tmp=np.transpose(tmp)
        ##For all snapshots average over the bound stars
        for j in range(len(tmp)):
            mean[i, j]=np.average(tmp[j][tmp[j]<1.])

    ##Averaging over all trials--this array should have a length equal to the number of snapshots.
    return np.average(mean, axis=0)

    test=get_eccs3(3,90,mypath2)
#New Plot
plt.plot(test)
plt.title('Control Plot')
plt.savefig('Control_Plot')
