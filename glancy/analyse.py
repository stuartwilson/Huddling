import h5py
import numpy as np
import pylab as pl

h5f = h5py.File('logs/agents.h5','r')

x = h5f['a1_x'][:]
y = h5f['a1_y'][:]

h5f.close()

F = pl.figure()
f = F.add_subplot(111)
f.plot(x,y)

pl.show()
