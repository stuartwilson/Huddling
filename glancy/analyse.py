import h5py
import numpy as np
import pylab as pl

h5f = h5py.File('logs/agents.h5','r')

foc = 0
x = h5f['a'+str(foc)+'_x'][:]
y = h5f['a'+str(foc)+'_y'][:]
b = h5f['a'+str(foc)+'_b'][:]

h5f.close()

F = pl.figure()
f = F.add_subplot(121)
f.plot(x,y)

f = F.add_subplot(122)
f.plot(b)

pl.show()
