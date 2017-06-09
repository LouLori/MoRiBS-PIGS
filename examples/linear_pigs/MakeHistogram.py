#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from os import system
from sys import argv

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
molecule            = "HF"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2
numbblocks	        = 20000                                                        #change param3
numbmolecules       = 1                                                            #change param4
tau                 = 0.001
dipolemoment        = 1.86                                                         #change param7
Rpt                 = 10.0

nrange              = 51		  						                           #change param11
skip                = 5

file1_name          = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
file1_name         += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"

file2_name          = ""                                                   #change param10
value_min           = 1                                                            #change param12

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13

nbins               = 100          #change


def plothistogram(srcfile, numbmolecules, molecule, dipolemoment, beta, tau, Rpt):
	col_costheta1, col_theta1 = loadtxt(srcfile, unpack=True, usecols=[0,1,2,3])

	hist1, bins1 = np.histogram(col_costheta1,nbins,density=True)
	print len(hist1), len(bins1)
	print np.sum(hist1)
	print np.sum(hist1*np.diff(bins1))

	plt.hist( col_costheta1, bins='auto', normed = 1)  # plt.hist passes it's arguments to np.histogram
	plt.grid(True)


	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$'+' R = '+str(Rpt)+r'$\AA$', fontsize =10)
	plt.xlabel(r'$\cos(\theta)$', fontsize = 20)
	plt.ylabel('Density', fontsize = 20)
	plt.xlim(-1.01,1.01)
	plt.grid(True)
	plt.show()

# Loop over jobs
for i in range(nrange):                                                  #change param13

	if (i>1 and i % skip == 0 ):

		if i % 2 != 0:
			value        = i
		else:
			value        = i+value_min

		numbbeads    = value
		beta         = tau*(numbbeads - 1)
		folder_run   = file1_name+str(value)+file2_name
		dest_dir     = dest_path + folder_run + "/results"
		src_file     = dest_dir+"/pigs_instant.eng"
		plothistogram(src_file, numbmolecules, molecule, dipolemoment, beta, tau, Rpt)
