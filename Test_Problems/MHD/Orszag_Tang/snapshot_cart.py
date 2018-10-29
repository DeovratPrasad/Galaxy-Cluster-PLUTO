import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
sys.path.append("/Users/prateek/Desktop/codes/Galaxy-Cluster-PLUTO/Tools/pyPLUTO/pyPLUTO/") #access linux environment variable
import pyPLUTO as pp

# script to make 2-D snapshots of various quantities

D=pp.pload(0)
x = D.x1
y = D.x2

for i in range(0,1):
	#plt.clf()
	D=pp.pload(i)
	plt.subplot(111)
	plt.contourf(x, y, D.vx1)
	#plt.xlim(0.0, 5.0)
	#plt.ylim(-5.0, 5.0)
	plt.axes().set_aspect('equal')
	plt.colorbar()
	plt.show()
	#fname = 'log10d_%03d.png'%i
	#print 'Saving frame', i
	#plt.savefig(fname)
