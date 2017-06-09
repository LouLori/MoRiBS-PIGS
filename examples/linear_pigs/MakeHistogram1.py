#!/usr/bin/env python
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
import time
from os import system
import os
import decimal
import support

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
molecule            = "HF"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2

Rpt1                = 7.0                                                         #change param6
dipolemoment        = 1.86                                                         #change param7

numbblocks          = 5000                                                        #change param3
numbmolecules       = 2                                                            #change param4

#tau                 = 0.002
beta                = 0.2

var1                = "beta"                                                       #change param10
var2                = "tau"

num1                = 2                                                            #change param12
numbbeads           = 31

#beta                = tau*(numbbeads-1)
tau                 = beta/(numbbeads-1)


ls1                 = '-'
ls2                 = '--'
color1              = 'blue'
color2              = 'red'


if  ((var1 == 'tau') or (var1 == 'beta')): 
	if var1 == "tau":
		value2 = beta
		value1 = tau
	else:
		value1 = beta
		value2 = tau
	Figfile      = "Figure-Density-fixed-"+var1+str(value1)+"Kinv-"
	Figfile     += var2+str(value2)+"Kinv-Rpt"+str(Rpt1)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	Figfile     += "-System"+str(numbmolecules)+str(molecule)+".png"


#======================================================================================================================
def plotenergy(Rpt, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1):

#
	#file1_name   = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name   = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name  += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"
	file2_name   = ""

	dest_path    = "/work/tapas/linear_rotors/"                                 
	folder_run   = file1_name+str(numbbeads)+file2_name
	dest_dir     = dest_path + folder_run
	srcfile      = dest_dir+"/results/pigs_instant.dof"

	dir_exact1   = '/home/tapas/CodesForEigenValues/DiagLinearRotorDimer/'
	dir_exact1  += 'Results-Jrot5-Rpt'+str(Rpt)+'Angstrom/'
#
	data         = loadtxt(srcfile,unpack=True, usecols=[numb_col])

	#plt.hist(data, bins='auto', normed = 1)  # plt.hist passes it's arguments to np.histogram
	if (Rpt == 5.0):
		plt.hist(data, bins=50, histtype='stepfilled', normed=True, color = color1, label = 'PIGS')
	else:
		plt.hist(data, bins=50, histtype='stepfilled', normed=True, color = color1, alpha = 0.5, label = 'PIGS')

	plt.grid(True)
	plt.ylabel('Density', fontsize = 20)
	if (numb_col == 3):
		plt.xlim(-.01,6.29)
		plt.xlabel(r'$\phi$ (radian)', fontsize = 20)
		src_exact = 'DensityPhi.txt'

		dir_exact   = dir_exact1 + src_exact
		col1, col2   = loadtxt(dir_exact,unpack=True, usecols=[0, 1])
		plt.plot(col1, col2, linestyle = ls1, color = 'black', lw = 3, label = 'Exact')
	else:
		plt.xlim(-1.01,1.01)

	if (numb_col == 0):
		plt.xlabel(r'$\mathrm{\vec{e}_{1} \cdot \vec{e}_{2} }$', fontsize = 20)
		src_exact = 'DensityCosThetaTilde.txt'

	if (numb_col == 1):
		plt.xlabel(r'$\cos(\theta_{1})$', fontsize = 20)
		src_exact = 'DensityCosTheta1.txt'

		dir_exact   = dir_exact1 + src_exact
		col1, col2   = loadtxt(dir_exact,unpack=True, usecols=[0, 1])
		plt.plot(col1, col2, linestyle = ls1, color = 'black', lw = 3, label = 'Exact')

	if (numb_col == 2):
		plt.xlabel(r'$\cos(\theta_{2})$', fontsize = 20)
		src_exact = 'DensityCosTheta2.txt'

		dir_exact   = dir_exact1 + src_exact
		col1, col2   = loadtxt(dir_exact,unpack=True, usecols=[0, 1])
		plt.plot(col1, col2, linestyle = ls1, color = 'black', lw = 3, label = 'Exact')
		

fig = plt.figure(figsize=(8, 4), dpi=100)

'''
if var1 == 'tau':
	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt1)+' '+r'$\AA$' )
if var1 == 'beta':
	plt.suptitle('Parameters: System '+str(numbmolecules)+' '+str(molecule)+', '+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt1)+' '+r'$\AA$' )

'''

#=========================================
#
# Fig1
#
#=========================================
num2   = 1
plt.subplot(num1, 2, num2)

numb_col = 0
plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1)
#plotenergy(Rpt2, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls2, color2)

#=========================================
#
# Fig2
#
#=========================================
num2   = 2
plt.subplot(num1, 2, num2)

numb_col = 1
plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1)
#plotenergy(Rpt2, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls2, color2)

#=========================================
#
# Fig2
#
#=========================================
num2   = 3
plt.subplot(num1, 2, num2)

numb_col = 2
plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1)
#plotenergy(Rpt2, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls2, color2)

#=========================================
#
# Fig4
#
#=========================================
num2   = 4
plt.subplot(num1, 2, num2)

numb_col = 3
plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1)
#plotenergy(Rpt2, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls2, color2)

#===============================================================================

plt.subplots_adjust(top=0.90, bottom=0.20, left=0.15, right=0.98, hspace=0.6,
                    wspace=0.5)
plt.legend(loc=3, bbox_to_anchor=(-0.80,-0.70), ncol=3, borderaxespad=0.)

fig.savefig(Figfile, dpi=100, format='png')
plt.show()
