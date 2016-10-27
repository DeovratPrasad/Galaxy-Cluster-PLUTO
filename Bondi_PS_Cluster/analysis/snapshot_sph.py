import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
sys.path.append("/home/prateek/Desktop/Public_Codes/Galaxy-Cluster-PLUTO/Tools/pyPLUTO/pyPLUTO")
import pyPLUTO as pp

# script to make 2-D snapshots of various quantities

D=pp.pload(0, w_dir="/home/prateek/Desktop/Public_Codes/Galaxy-Cluster-PLUTO/Bondi_PS_Cluster/")
x = np.outer(D.x1,np.sin(D.x2))
z = np.outer(D.x1,np.cos(D.x2))

for i in range(0,10):
	plt.clf()
	D=pp.pload(i, w_dir="/home/prateek/Desktop/Public_Codes/Galaxy-Cluster-PLUTO/Bondi_PS_Cluster/")
	plt.subplot(111)
	plt.contourf(x, z, np.log10(D.rho), np.linspace(-4.0,0.0, 40))
	plt.xlim(0.0, 5.0)
	plt.ylim(-5.0, 5.0)
	plt.axes().set_aspect('equal')
	plt.colorbar()
	plt.show()
	fname = 'log10d_%03d.png'%i
	print 'Saving frame', i
	plt.savefig(fname)
