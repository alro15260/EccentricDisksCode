%pylab inline
import ipyvolume as ipv
import ipyvolume.pylab as p3
import numpy as np
import rebound
from rebound import hash as h

from astropy.io import ascii

p3.clear()
##Location of npz files produced by other script
loc="./"
x=np.load(loc+'x.npz')
y=np.load(loc+'y.npz')
z=np.load(loc+'z.npz')

s2=p3.plot(x['arr_0'], y['arr_0'], z['arr_0'])
# s3=p3.plot(x['arr_0'][:, 1, :], y['arr_0'][:, 1, :], z['arr_0'][:, 1, :], color='red', size=0.3, alpha='0.5')

ipv.pylab.xlim(-1,1)
ipv.pylab.ylim(-1,1)
ipv.pylab.zlim(-1,1)
ipv.pylab.view(azimuth=None, elevation=-90, distance=1.5)

p3.show()
ipv.animation_control(s2)
