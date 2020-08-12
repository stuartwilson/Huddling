import h5py
import numpy as np
import pylab as pl

h5f = h5py.File('test/agents.h5','r')

x1 = h5f['a1_x'][:]
y1 = h5f['a1_y'][:]

x2 = h5f['a2_x'][:]
y2 = h5f['a2_y'][:]


h5f.close()

F = pl.figure()
f = F.add_subplot(111)
f.plot(x1,y1)
f.plot(x2,y2)

pl.show()
