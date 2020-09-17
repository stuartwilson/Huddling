import h5py
import numpy as np
import pylab as pl

h5f = h5py.File('test/agents.h5','r')

i = 11
x = h5f['a'+str(i)+'_x'][:]
y = h5f['a'+str(i)+'_y'][:]
w1 = h5f['a'+str(i)+'_w_1'][:]
w2 = h5f['a'+str(i)+'_w_0'][:]

h5f.close()

F = pl.figure(figsize=(12,6))
f = F.add_subplot(121)
f.plot(x[-1000:],y[-1000:])
f.axis(np.array([-1,1,-1,1.])*10)
f = F.add_subplot(122)
f.plot(w1)
f.plot(w2)

pl.show()
