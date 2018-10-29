import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.environ['PLUTO_DIR']+"/Tools/pyPLUTO/pyPLUTO") #access linux environment variable
import pyPLUTO as pp
import itertools

# script to make 1-D snapshots of n_e, TkeV, emission/volume/mass weighted

def cf_TN(TkeV):
	lam = 1.e-22*( 8.6e-3*TkeV**-1.7 + 0.058*math.sqrt(TkeV) + 0.063 )
	return lam 

mu = 0.62; mue = 1.17;  unit_d = 1.67262171e-24; unit_l = 3.0856775807e18*1.e5; unit_v = 1.e8
eV = 1.16e4; gamm = 5./3.; gamm1 = gamm-1.0; kB = 1.3806505e-16; mp = 1.67262171e-24

mui = 1./(1./mu-1./mue)

option = 'm' #valid option e/v/m

D=pp.pload(0, w_dir="/home/prateek/Desktop/Public_Codes/Galaxy-Cluster-PLUTO/Bondi_PS_Cluster/")
r = D.x1; th = D.x2; dr = D.dx1; dth = D.dx2; dph = D.dx3
nr = np.size(D.x1); nt = np.size(D.x2); np = np.size(D.x3)

for l in range(0,10):
	fname = open('d_vs_r_'+option+str(l)+'.dat','w')
        D=pp.pload(l, w_dir="/home/prateek/Desktop/Public_Codes/Galaxy-Cluster-PLUTO/Bondi_PS_Cluster/")
	for i in range(nr):
		den_ne = 0.0; vol = 0.0; Temp_keV = 0.0
		for j, k in itertools.product(range(nt), range(np)):
			dvol = r[i]*r[i]*dr[i]*math.sin(th[j])*dth[j]*dph[k]
			TkeV = D.prs[i,j]*unit_v*unit_v*mu*mp/(D.rho[i,j]*kB*eV*1.e3)
			n_e = D.rho[i,j]*unit_d/(mue*mp); n_i = D.rho[i,j]*unit_d/(mui*mp)
			if (option == "e"):
				if ( TkeV > 0.1 ):
					dvol = dvol*n_e*n_i*cf_TN(TkeV)
					vol = vol + dvol
					den_ne = den_ne + n_e*dvol
					Temp_keV = Temp_keV + TkeV*dvol
			if (option == "v"):
				vol = vol + dvol
				den_ne = den_ne + n_e*dvol
				Temp_keV = Temp_keV + TkeV*dvol
			if (option == "m"):
				dvol = dvol*D.rho[i,j]
                                vol = vol + dvol
                                den_ne = den_ne + n_e*dvol
                                Temp_keV = Temp_keV + TkeV*dvol
			
		den_ne = den_ne/vol; Temp_keV = Temp_keV/vol
		fname.write('%.10e %.10e %.10e\n'% (r[i], den_ne, Temp_keV))
	fname.close()	
	print 'Saving frame', l
