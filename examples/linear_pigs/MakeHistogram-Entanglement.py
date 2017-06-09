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

Rpt2                = 5.0                                                         #change param6
Rpt1                = 10.0                                                         #change param6
dipolemoment        = 1.86                                                         #change param7

numbblocks          = 5000                                                        #change param3
numbmolecules       = 4                                                            #change param4

tau                 = 0.002

var1                = "beta"
var2                = "tau"                                                       #change param10

num1                = 2                                                            #change param12
numbbeads           = 11

beta                = tau*(numbbeads-1)


ls1                 = '-'
ls2                 = '--'
color1              = 'blue'
color2              = 'red'


#======================================================================================================================
def plotenergy(Rpt, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell):

#
	file1_name   = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name  += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"
	file2_name   = ""

	dest_path    = "/work/tapas/linear_rotors/"                                 
	folder_run   = file1_name+str(numbbeads)+file2_name
	dest_dir     = dest_path + folder_run
	srcfile      = dest_dir+"/results/pigs_instant.dof"

	data         = loadtxt(srcfile,unpack=True, usecols=[numb_col])

	plt.hist(data, bins=50, histtype='stepfilled', normed=True, color = color1, alpha = 0.5, label = 'PIGS')

	plt.grid(True)
	plt.ylabel('Density', fontsize = 20)
	plt.xlabel(xlabell, fontsize = 20)

for i in range(2):

	if (i == 0):
		bead = "(M - 1) Bead"
		beadw = "M1"
	if (i == 1):
		bead = "M Bead"
		beadw = "M"

	for j in range(3):

		fig = plt.figure(figsize=(8, 4), dpi=100)

		if var1 == 'tau':
			plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt1)+' '+r'$\AA$' )
		if var1 == 'beta':
			plt.suptitle('Parameters: System '+str(numbmolecules)+' '+str(molecule)+', '+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt1)+' '+r'$\AA$  '+bead)

		if  ((var1 == 'tau') or (var1 == 'beta')): 
			if var1 == "tau":
				value2 = beta
				value1 = tau
			else:
				value1 = beta
				value2 = tau

		if (j == 0):
			cord = "X"
		if (j == 1):
			cord = "Y"
		if (j == 2):
			cord = "Z"


		Figfile      = "Figure-Entanglement-pot0-Density-Of-"+cord+"-Bead-"+beadw+"-fixed-"+var1+str(value1)+"Kinv-"
		Figfile     += var2+str(value2)+"Kinv-Rpt"+str(Rpt1)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		Figfile     += "-System"+str(numbmolecules)+str(molecule)+".png"

		if (i == 0):
			jj = 3*i + j*4
		else:
			jj = 3*i + j*4 + 9
			
#=========================================
#
# Fig1
#
#=========================================
		num2   = 1
		plt.subplot(num1, 2, num2)

		numb_col = jj
		subscript = 1
		xlabell = cord+str(subscript)
		plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell)

#=========================================
#
# Fig2
#
#=========================================
		num2   = 2
		plt.subplot(num1, 2, num2)

		numb_col = jj + 1
		subscript = 2
		xlabell = cord+str(subscript)
		plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell)

#=========================================
#
# Fig2
#
#=========================================
		num2   = 3
		plt.subplot(num1, 2, num2)

		numb_col = jj + 2
		subscript = 3
		xlabell = cord+str(subscript)
		plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell)

#=========================================
#
# Fig4
#
#=========================================
		num2   = 4
		plt.subplot(num1, 2, num2)

		numb_col = jj + 3
		subscript = 4
		xlabell = cord+str(subscript)
		plotenergy(Rpt1, dipolemoment, tau, numbblocks, numbmolecules, molecule, numbbeads, numb_col, ls1, color1, xlabell)

#===============================================================================

		plt.subplots_adjust(top=0.90, bottom=0.20, left=0.15, right=0.98, hspace=0.6, wspace=0.5)
		plt.legend(loc=3, bbox_to_anchor=(-0.80,-0.70), ncol=3, borderaxespad=0.)

		fig.savefig(Figfile, dpi=100, format='png')
#plt.show()
