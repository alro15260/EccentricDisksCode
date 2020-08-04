from astropy.table import Table
import matplotlib.pyplot as plt
import sys
import numpy as np

from scipy.optimize import curve_fit

def cum_dist(x, xmin, xmax, p):
    return (x**(1-p)-xmin**(1-p))/(xmax**(1-p)-xmin**(1-p))


elem='a'
names=np.genfromtxt('names', dtype=str)

f=open('pfit', 'w')
for ii in range(0, 4994, 10):
	elem_list=np.concatenate([Table.read(name.replace('.bin', '_elems.hdf5'), '/{0}'.format(ii))[elem] for name in names])
	elem_list=elem_list[elem_list>0]
	np.savetxt('a_list_{0}.txt'.format(ii), elem_list)

	xx=np.sort(elem_list)
	yy=(np.array(range(len(xx)))).astype(float)/float(len(xx))
	popt,pcov=curve_fit(cum_dist, xx, yy, [1., 2., 1.1])
	f.write('{0} {1} {2} {3}\n'.format(ii, popt[0], popt[1], popt[2]))

f.close()
