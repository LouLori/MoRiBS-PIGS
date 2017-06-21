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
#TypeCal             = 'PIGS'
TypeCal             = 'ENT'

molecule            = "HF"                                                     
molecule_rot        = "HF"                                                   

numbblocks	        = 40000
numbmolecules       = 2
numbpass            = 10

tau                 = 0.001                                               

Rpt                 = 10.05
dipolemoment        = 1.86

particleA           = 1

#ENT_TYPE = "SWAP"
ENT_TYPE = "BROKENPATH"
#ENT_TYPE = "REGULARPATH"

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-e0vsbeads-"+ENT_TYPE 


file2_name          = ""                                                           #change param13

# 

var1                = "beta"
var2                = "tau"                                                       #change param10


numbbeads           = 5
beta                = tau*(numbbeads-1)


ls1                 = '-'
ls2                 = '--'
color1              = 'blue'
color2              = 'red'


#======================================================================================================================
def plotenergy(Rpt, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell):

#
	dest_path    = "/work/tapas/linear_rotors/"                                 
	folder_run   = file1_name+str(numbbeads)+file2_name
	dest_dir     = dest_path + folder_run
	srcfile      = dest_dir+"/results/pigs_instant.dof"

	data         = loadtxt(srcfile,unpack=True, usecols=[numb_col])

#	plt.hist(data, bins=50, histtype='stepfilled', normed=True, color = color1, alpha = 0.5, label = 'PIGS')
	plt.hist(data, bins=50, histtype='stepfilled', normed=True, color = color1, label = 'PIGS')

	plt.grid(True)
	plt.ylabel('Density', fontsize = 20)
	plt.xlabel(xlabell, fontsize = 20)

fig = plt.figure(figsize=(8, 4), dpi=100)

if var1 == 'tau':
	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt)+' '+r'$\AA$' )
if var1 == 'beta':
	plt.suptitle('Parameters: System '+str(numbmolecules)+' '+str(molecule)+', '+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt)+' '+r'$\AA$  ')

if  ((var1 == 'tau') or (var1 == 'beta')): 
	if var1 == "tau":
		value2 = beta
		value1 = tau
	else:
		value1 = beta
		value2 = tau

cord = "Estimator"

if ((ENT_TYPE == 'SWAP') or (ENT_TYPE == 'REGULARPATH')):
	jj = 35 
	num1 = 1                                                            #change param12

if (ENT_TYPE == 'SWAP'):
	xlabel1 = "DM/NM"

if (ENT_TYPE == 'REGULARPATH'):
	xlabel1 = "NM/DM"

if (ENT_TYPE == 'BROKENPATH'):
	jj = 33 
	num1 = 2                                                            #change param12
	xlabel1 = "NM"
	xlabel2 = "DM"


Figfile      = "Figure-Entanglement-pot0-Density-Of-"+cord+"-fixed-"+var1+str(value1)+"Kinv-"
Figfile     += var2+str(value2)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
Figfile     += "-System"+str(numbmolecules)+str(molecule)+".png"

			
#=========================================
#
# Fig1
#
#=========================================
num2   = 1
plt.subplot(num1, 1, num2)

numb_col = jj
plotenergy(Rpt, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabel1)

#=========================================
#
# Fig2
#
#=========================================
if (ENT_TYPE == 'BROKENPATH'):
	num2   = 2
	plt.subplot(num1, 1, num2)

	numb_col = jj + 1
	plotenergy(Rpt, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabel2)

#===============================================================================

plt.subplots_adjust(top=0.90, bottom=0.20, left=0.15, right=0.98, hspace=0.6, wspace=0.5)
plt.legend(loc=3, bbox_to_anchor=(-0.80,-0.70), ncol=3, borderaxespad=0.)

fig.savefig(Figfile, dpi=100, format='png')
plt.show()
