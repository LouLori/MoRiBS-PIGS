#!/usr/bin/env python
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
#molecule            = "HF-C60"                                                     #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                         #change param2
molecule_rot        = "HF"                                                         #change param2

Rpt                 = 7.0                                                         #change param6
dipolemoment        = 1.86                                                         #change param7
nrange              = 51                                                            #change param5

numbblocks          = 40000                                                        #change param3
numbmolecules       = 2                                                            #change param4
numbbeads           = 129

beta                = 0.2                                                        #change param7
tau                 = 0.002

Avg_total_energy    = -11.823                                                      #change param8
Avg_potential_energy= -21.1273                                                     #change param9
Avg_total_energy    = -0.0115404                                                   #change param8
Avg_potential_energy= -0.0230796                                                   #change param9
Avg_rotational_energy = 0.0115393
Avg_costheta        = -0.0159729

var2                = "beta"  #varaible                                                     #change param10
var1                = "tau" # fixed
var3                = "tau"

num1                = 2                                                            #change param12
trunc 				= 40000                                                            #change param13
trunc1              = 15
font=18

if Rpt == 10.0:
	if ((var1 == 'tau') or (var1 == 'beta') and (Rpt == 10.0) and (dipolemoment == 1.86)):
		if numbmolecules == 1:
			Avg_total_energy      = -12.14311781
			Avg_potential_energy  = -21.64771613
			Avg_rotational_energy = 9.504598325
			Avg_costheta          = 0.4319538676
			xmin = 1.9
			xmax = 10.1
			yminTot = -6000.0
			ymaxTot = 0.0
			yminPot = yminTot
			ymaxPot = ymaxTot
			yminRot = -1.0
			ymaxRot = 300.0
			yminAng = 0.4
			ymaxAng = 1.0
		if numbmolecules == 2:
			Avg_total_energy      = -3.424708242
			Avg_potential_energy  = -6.760526705
			Avg_rotational_energy = 3.335818463
			Avg_costheta          = 0.003822145993
			xmin = 1.9
			xmax = 10.1
			yminTot = -6000.0
			ymaxTot = 0.0
			yminPot = yminTot
			ymaxPot = ymaxTot
			yminRot = -1.0
			ymaxRot = 30.0
			yminAng = 0.0
			ymaxAng = 1.0
if  ((var1 == 'tau') or (var1 == 'beta')): 
	if var1 == "tau":
		value = beta
	else:
		value = tau
	Figfile      = "Figure-Energy-vs-"+var1+"-fixed-"
	Figfile     += var2+str(value)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	Figfile     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".png"


if Rpt == 7.0:
	if ((var1 == 'tau') or (var1 == 'beta') and (numbmolecules == 2) and (Rpt == 7.0) and (dipolemoment == 1.86)):
		Avg_total_energy      = -26.85431770
		Avg_potential_energy  = -49.57886525
		Avg_rotational_energy = 22.72454755   
		Avg_costheta          = 0.06758428991
		xmin = 1.9
		xmax = 10.1
		yminTot = -6000.0
		ymaxTot = 0.0
		yminPot = yminTot
		ymaxPot = ymaxTot
		yminRot = -1.0
		ymaxRot = 300.0
		yminAng = 0.4
		ymaxAng = 1.0

if ((var1 == 'Rpt') and (dipolemoment == 1.86)):
	tau     = beta/(numbbeads - 1)

	if numbmolecules == 1:
		xmin    = 1.9
		xmax    = 10.1
		yminTot = -6000.0
		ymaxTot = 0.0
		yminPot = yminTot
		ymaxPot = ymaxTot
		yminRot = -1.0
		ymaxRot = 300.0
		yminAng = 0.4
		ymaxAng = 1.0
	if numbmolecules == 2:
		xmin    = 1.9
		xmax    = 10.1
		yminTot = -6000.0
		ymaxTot = 0.0
		yminPot = yminTot
		ymaxPot = ymaxTot
		yminRot = -1.0
		ymaxRot = 300.0
		yminAng = 0.0
		ymaxAng = 1.0


#======================================================================================================================

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('--','-','-.',':'))
		self.color2               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('--','-','-.',':'))

def plotenergy(status,func,var1,var2,var3,Rpt,dipolemoment,beta,tau,numbmolecules,molecule,numbblocks,x,trunc,trunc1,numbbeads):
	markersize_fig = 1

	if ((var1 == "tau") or (var1 == "beta")):
		if var1 == "tau":
			value = beta
		else:
			value = tau
		if status == "Rotation":
			srcfile      = "AngularDOF-vs-"+var1+"-fixed-"
			srcfile     += var2+str(value)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
			srcfile     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

			data         = loadtxt(srcfile,unpack=True, usecols=[1,2,3,4,5])
			var, costheta, theta, err_costheta, err_theta = data
			var          = var[0:trunc1]
			costheta     = costheta[0:trunc1]
			theta        = theta[0:trunc1]
			err_costheta = err_costheta[0:trunc1]
			err_theta    = err_theta[0:trunc1]

		if status == "Energy":
			srcfile      = "Energy-vs-"+var1+"-fixed-"
			srcfile     += var2+str(value)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
			srcfile     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

			data         = loadtxt(srcfile,unpack=True, usecols=[1,2,3,4,5,6,7])
			var, pot, tot, rot, err_pot, err_tot, err_rot = data
			var          = var[0:trunc1]
			pot          = pot[0:trunc1]
			tot          = tot[0:trunc1]
			rot          = rot[0:trunc1]
			err_pot      = err_pot[0:trunc1]
			err_tot      = err_tot[0:trunc1]
			err_rot      = err_rot[0:trunc1]

	if var1 == "Rpt":
		srcexact         = "EigenValuesFor"+str(numbmolecules)+molecule+"-DipoleMoment"+str(dipolemoment)+".txt"
		dataexact        = loadtxt(srcexact,unpack=True, usecols=[0,1,3,4,5])
		exactRpt, exactTot, exactPot, exactRot, exactAng = dataexact

		if status == "Rotation":
			srcfile      = "AngularDOF-vs-"+var1+"-fixed-"
			srcfile     += var2+str(beta)+"Kinv-"+var3+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
			srcfile     += "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-trunc"+str(trunc)+".txt"

			data         = loadtxt(srcfile,unpack=True, usecols=[1,2,3,4,5])
			var, costheta, theta, err_costheta, err_theta = data
			var          = var[0:trunc1]
			costheta     = costheta[0:trunc1]
			theta        = theta[0:trunc1]
			err_costheta = err_costheta[0:trunc1]
			err_theta    = err_theta[0:trunc1]

		if status == "Energy":
			srcfile      = "Energy-vs-"+var1+"-fixed-"
			srcfile     += var2+str(beta)+"Kinv-"+var3+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
			srcfile     += "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-trunc"+str(trunc)+".txt"

			data         = loadtxt(srcfile,unpack=True, usecols=[1,2,3,4,5,6,7])
			var, pot, tot, rot, err_pot, err_tot, err_rot = data
			var          = var[0:trunc1]
			pot          = pot[0:trunc1]
			tot          = tot[0:trunc1]
			rot          = rot[0:trunc1]
			err_pot      = err_pot[0:trunc1]
			err_tot      = err_tot[0:trunc1]
			err_rot      = err_rot[0:trunc1]

	if func == 'TotalEnergy':
		val     = tot
		err_val = err_tot

		if var1 == "Rpt":
			plt.plot(exactRpt, exactTot, linestyle = '-', color = 'black', label = 'Exact', lw = 3)
			
	if func == 'PotentialEnergy':
		val     = pot
		err_val = err_pot

		if var1 == "Rpt":
			plt.plot(exactRpt, exactPot, linestyle = '-', color = 'black', label = 'Exact', lw = 3)

	if func == 'RotationalEnergy':
		val     = rot
		err_val = err_rot

		if var1 == "Rpt":
			plt.plot(exactRpt, exactRot, linestyle = '-', color = 'black', label = 'Exact', lw = 3)

	if func == "RelativeAngle":
		val     = costheta
		err_val = err_costheta

		if var1 == "Rpt":
			plt.plot(exactRpt, exactAng, linestyle = '-', color = 'black', label = 'Exact', lw = 3)

		
	plt.plot(var, val, linestyle = x.line1.next(), color = x.color1.next(), marker = x.marker1.next(), label = 'PIGS', lw = 1)
	plt.errorbar(var, val, yerr = err_val, linestyle="None", marker="None", color= x.color2.next())

def fitFunc(var, a, b):
	return a + b*var*var

def plotfitting(status, var1, var2, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc, x, trunc1,  b):
	markersize_fig = 10
	if status == "Energy":
		srcfile      = "Energy-vs-"+var1+"-fixed-"
		srcfile     += var2+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		srcfile     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		data        = loadtxt(srcfile, unpack=True, usecols=[1,2,3,5,6])
		var, pot, tot, err_pot, err_tot = data
		var = var[0:trunc1]
		pot = pot[0:trunc1]
		tot = tot[0:trunc1]
		err_pot = err_pot[0:trunc1]
		err_tot = err_tot[0:trunc1]

	if func == 'TotalEnergy':
		val     = tot
		err_val = err_tot
		fitParams, fitCovariances = curve_fit(fitFunc, var, val, sigma = err_val)
		plt.plot(var, fitFunc(var, fitParams[0], fitParams[1]), linestyle = x.line1.next(), color = x.color1.next(), marker = x.marker1.next(), label = 'Fit', lw = 1)
	if func == 'PotentialEnergy':
		val     = pot
		err_val = err_pot
		fitParams, fitCovariances = curve_fit(fitFunc, var, val, sigma = err_val)
		plt.plot(var, fitFunc(var, fitParams[0], fitParams[1]), linestyle = x.line1.next(), color = x.color1.next(), marker = x.marker1.next(), label = 'Fit', lw = 1)
		
#	fitParams, fitCovariances = curve_fit(fitFunc, var, val, sigma = err_val)
#	plt.plot(var, fitFunc(var, fitParams[0], fitParams[1]), linestyle = x.line1.next(), color = x.color1.next(), marker = x.marker1.next(), label = 'Fit', lw = 1)


def plotenergylabel(num1,num2,Avg_energy,xmin,xmax,ymin,ymax,xlabel1,ylabel1,font):
	if xlabel1 != "Rpt":
		plt.axhline(y=Avg_energy, color='black', lw = 3.0, linestyle='--', label = 'Exact')

	plt.grid(True)
#	plt.xlim(0,0.2)
#	plt.legend(loc='upper right')


	if xlabel1 == "Rpt":
		plt.xlabel(r'R ($\AA$)', fontsize = font)
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
	if xlabel1 == "tau":
		plt.xlabel(r'$\mathrm{\tau (K^{-1})}$', fontsize = font)
	if xlabel1 == "beta":
		plt.xlabel(r'$\mathrm{\beta (K^{-1})}$', fontsize = font)

	if ylabel1 == "TotalEnergy":
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (Kelvin)}$', fontsize = font)
	if ylabel1 == "PotentialEnergy":
		plt.ylabel(r'$\mathrm{\langle V_{0} \rangle (Kelvin)}$', fontsize = font)
	if ylabel1 == "RotationalEnergy":
		plt.ylabel(r'$\mathrm{\langle T_{0} \rangle (Kelvin)}$', fontsize = font)
		#plt.ylim(-1,20)
	if ylabel1 == "RelativeAngle":
		plt.ylabel(r'$\mathrm{\langle \vec{e}_{1} \cdot \vec{e}_{2} \rangle}$', fontsize = font)

	#plt.xticks(np.linspace(-0, 0.2, 5, endpoint=True))
	#plt.text(0.1, -2.0, r'$\tau$ = 0.001 K $^{-1}$')


fig = plt.figure(figsize=(8, 4), dpi=100)

'''
if var1 == 'Rpt':
	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$')
if var1 == 'tau':
	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, Rpt = '+str(Rpt)+' '+r'$\AA$')
if var1 == 'beta':
	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$, Rpt = '+str(Rpt)+' '+r'$\AA$')
'''

#=========================================
#
# Fig1
#
#=========================================
num2   = 1
x      = Specification()
status = "Energy"
func   = "TotalEnergy"
plt.subplot(num1, 2, num2)

plotenergy(status,func,var1,var2,var3,Rpt,dipolemoment,beta,tau,numbmolecules,molecule,numbblocks,x,trunc,trunc1,numbbeads)

# fitting done here
if var1 == "tau":
	b = -100.0
	plotfitting(status, var1, var2, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc, x, trunc1, b)

xlabel = var1
ylabel = func

plotenergylabel(num1,num2,Avg_total_energy,xmin,xmax,yminTot,ymaxTot,xlabel,ylabel,font)

#=========================================
#
# Fig2
#
#=========================================
num2   = 2
x      = Specification()
status = "Energy"
func   = "PotentialEnergy"
plt.subplot(num1, 2, num2)

plotenergy(status,func,var1,var2,var3,Rpt,dipolemoment,beta,tau,numbmolecules,molecule,numbblocks,x,trunc,trunc1,numbbeads)

# fitting done here
if var1 == "tau":
	b = -100.0
	plotfitting(status, var1, var2, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc, x, trunc1, b)

xlabel = var1
ylabel = func
plotenergylabel(num1,num2,Avg_potential_energy,xmin,xmax,yminPot,ymaxPot,xlabel,ylabel,font)

plt.legend(loc=2, bbox_to_anchor=(-1.05,-2.15), ncol=3, borderaxespad=0.)
#=========================================
#
# Fig3
#
#=========================================
num2   = 3
x      = Specification()
status = "Energy"
func   = "RotationalEnergy"
plt.subplot(num1, 2, num2)

plotenergy(status,func,var1,var2,var3,Rpt,dipolemoment,beta,tau,numbmolecules,molecule,numbblocks,x,trunc,trunc1,numbbeads)

xlabel = var1
ylabel = func

plotenergylabel(num1,num2,Avg_rotational_energy,xmin,xmax,yminRot,ymaxRot,xlabel,ylabel,font)

#=========================================
#
# Fig3
#
#=========================================
num2   = 4
x      = Specification()
status = "Rotation"
func   = "RelativeAngle"
plt.subplot(num1, 2, num2)

plotenergy(status,func,var1,var2,var3,Rpt,dipolemoment,beta,tau,numbmolecules,molecule,numbblocks,x,trunc,trunc1,numbbeads)

xlabel = var1
ylabel = func

plotenergylabel(num1,num2,Avg_costheta,xmin,xmax,yminAng,ymaxAng,xlabel,ylabel,font)

#===============================================================================

plt.subplots_adjust(top=0.95, bottom=0.25, left=0.15, right=0.98, hspace=0.6,
                    wspace=0.5)
#plt.legend(loc=3, bbox_to_anchor=(-0.75,-0.70), ncol=5, borderaxespad=0.)

fig.savefig(Figfile, dpi=100, format='png')
#call(["convert", 'image_output.png', 'image_output.pdf'])
plt.show()
