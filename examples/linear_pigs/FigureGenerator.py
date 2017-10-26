#!/usr/bin/env python
import numpy as np
from numpy import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import support
import os

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('--','-','-.',':'))
		self.color2               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('--','-','-.',':'))

def	FigureENT(FileToBePlot,FilePlot,TypeCal,variableName,parameter,ExactValue,numbmolecules,molecule,Rpt,dipolemoment):
	outfile = FilePlot
	fig = plt.figure(figsize=(16, 8), dpi=400)

	font=22
	#plt.grid(True)

	#x      = Specification()

	beta1, purity1, entropy1, err_purity1, err_entropy1 = loadtxt(FileToBePlot,unpack=True, usecols=[1, 4, 5, 8, 9])
	plt.subplot(1,2, 1)
	plt.errorbar(beta1, purity1, yerr=err_purity1, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
	if (variableName == "beta"):
		data2        = loadtxt(ExactValue,unpack=True, usecols=[0,2])
		beta2, entropy2 = data2
		plt.xlim(0,0.201)
		plt.plot(beta2, entropy2, color = 'black', ls = '--', label = 'Diagonalization', linewidth=1)

	if (variableName == "tau"):
		plt.xlim(0,)
		#plt.axhline(y=ExactValue, color='black', lw = 2.0, linestyle='--', label = 'Diagonalization')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	print(textpositionx)
	print(textpositiony)


	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  (Kelvin^{-1}})$', fontsize = font)
	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)

	plt.ylabel(r'$\mathrm{\langle P_{A} \rangle}$', fontsize = font)
	plt.subplots_adjust(top=0.95, bottom=0.30, left=0.25, right=1.0, hspace=0.0, wspace=1.0)

	plt.subplot(1, 2, 2)
	error1 = err_entropy1
	plt.errorbar(beta1, entropy1, yerr=error1, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
	if (variableName == "beta"):
		data2        = loadtxt(ExactValue,unpack=True, usecols=[0,2])
		beta2, entropy2 = data2
		plt.xlim(0,0.201)
		plt.plot(beta2, entropy2, color = 'black', ls = '--', label = 'Diagonalization', linewidth=1)

	if (variableName == "tau"):
		plt.xlim(0,)
		plt.axhline(y=ExactValue, color='black', lw = 2.0, linestyle='--', label = 'Diagonalization')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	print(textpositionx)
	print(textpositiony)


	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  (Kelvin^{-1}})$', fontsize = font)
	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)
	
	plt.ylabel(r'$\mathrm{S_{2}}$', fontsize = font)


	plt.subplots_adjust(top=0.95, bottom=0.10, left=0.15, right=0.90, hspace=0.0, wspace=0.2)
	plt.legend(bbox_to_anchor=(0.30, 0.40), loc=2, borderaxespad=0., shadow=True, fontsize = font)
	plt.savefig(outfile, dpi = 200, format = 'pdf')

	call(["open", outfile])
	#call(["okular", outfile])
	#plt.show()
	'''
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
	'''
	'''
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
	'''

def	FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint):
	outfile = FilePlot
	fig = plt.figure(figsize=(6, 4), dpi=200)

	font=12
	datacorr = genfromtxt(FileToBePlot)
	FuncCorr = np.zeros((numbmolecules,numbmolecules))
	ErrorFuncCorr = np.zeros((numbmolecules,numbmolecules))
	ii = 0
	for i in range(numbmolecules-1):
		for j in range((i+1),numbmolecules):
			nc = 2+(2*ii)
			nec = nc+1 
			FuncCorr[i,j] = datacorr[nc]
			FuncCorr[j,i] = FuncCorr[i,j]
			ErrorFuncCorr[i,j] = datacorr[nec]
			ErrorFuncCorr[j,i] = ErrorFuncCorr[i,j]
			ii = ii+1
	val1 = np.arange(numbmolecules)
	val2 = np.zeros(numbmolecules)
	val3 = np.zeros(numbmolecules)
	for j in range(numbmolecules):
		if (j != RefPoint):
			val2[j] = FuncCorr[RefPoint:(RefPoint+1),j]
			val3[j] = ErrorFuncCorr[RefPoint:(RefPoint+1),j]
	
	plt.errorbar(val1[RefPoint+1:], val2[RefPoint+1:], yerr=val3[RefPoint+1:], color = 'b', ls = '-', label = 'PIGS', linewidth=1)
	print(val1[RefPoint+1:])

	plt.grid(True)
	#plt.xlim(1.98, 6.02)

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	plt.xlim(xmin-0.1,xmax+0.1)
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{j}$', fontsize = font)
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
		plt.text(textpositionx, textpositiony-5*deltay, r'$\mathrm{i} = $' + str(RefPoint), fontsize=10)

	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)
		plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony+2*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
		plt.text(textpositionx, textpositiony+0*deltay, r'$\mathrm{i} = $' + str(RefPoint), fontsize=10)
	
	if (TypeCorr == "Total"):
		plt.ylabel(r'$\mathrm{C_{ij}}$', fontsize = font)
	if (TypeCorr == "XCorr"):
		plt.ylabel(r'$\mathrm{C^{X}_{ij}}$', fontsize = font)
	if (TypeCorr == "YCorr"):
		plt.ylabel(r'$\mathrm{C^{Y}_{ij}}$', fontsize = font)
	if (TypeCorr == "ZCorr"):
		plt.ylabel(r'$\mathrm{C^{Z}_{ij}}$', fontsize = font)
	if (TypeCorr == "XYCorr"):
		plt.ylabel(r'$\mathrm{C^{X,Y}_{ij}}$', fontsize = font)

	plt.subplots_adjust(top=0.95, bottom=0.20, left=0.17, right=0.95, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.70, 0.70), loc=2, borderaxespad=0., shadow=True, fontsize = font)
	plt.savefig(outfile, dpi = 200, format = 'pdf')

	#call(["okular", outfile])
	call(["open", outfile])
	#plt.show()

def	FigureEnergyPIGS(FileToBePlot,FilePlot,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment):
	outfile = FilePlot
	fig = plt.figure(figsize=(4, 8), dpi=200)

	font=12
	valTau, valRotEnergy, valPotEnergy, valTotalEnergy, errorRotEnergy, errorPotEnergy, errorTotalEnergy = genfromtxt(FileToBePlot, unpack=True, usecols=[1, 3, 4, 5, 7, 8, 9])

	plt.subplot(3, 1, 1)
	YLabel = "Rotational"
	PlotEnergyPIGS(font, valTau,valRotEnergy,errorRotEnergy,variableName,YLabel)
	ymin, ymax = plt.ylim()
	xmin, xmax = plt.xlim()
	PlotLabelEnergyPIGS(font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)

	plt.subplot(3, 1, 2)
	YLabel = "Potential"
	PlotEnergyPIGS(font, valTau,valPotEnergy,errorPotEnergy,variableName,YLabel)

	plt.subplot(3, 1, 3)
	YLabel = "Total"
	PlotEnergyPIGS(font, valTau,valTotalEnergy,errorTotalEnergy,variableName,YLabel)

	plt.subplots_adjust(top=0.95, bottom=0.10, left=0.25, right=0.95, hspace=0.30, wspace=1.0)
	#plt.legend(bbox_to_anchor=(0.70, 0.70), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi = 200, format = 'pdf')

	call(["open", outfile])
	#call(["okular", outfile])
	#plt.show()

def PlotLabelEnergyPIGS(font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment):
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if (variableName == "beta"):
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)

	if (variableName == "tau"):
		plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony+2*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
	

def PlotEnergyPIGS(font,val1,val2,val3,variableName,YLabel):
	#plt.grid(True)
	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)
	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  (Kelvin^{-1}})$', fontsize = font)
		
	if (YLabel == "Total"):
		plt.ylabel(r'$\mathrm{E_{0}}$ (Kelvin)', fontsize = font)
		if (variableName == "tau"):
			plotfitting(val1, val2, val3)
			plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
		if (variableName == "beta"):
			plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
			ExactValueFile = "ResultsOfPIGSENT/ENT-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-Entropy-vs-beta-fixed-tau0.005Kinv-System2HF-ParticleA1-by-Dmitri.txt"
			data2        = loadtxt(ExactValueFile,unpack=True, usecols=[0,1])
			beta2, energy2 = data2
			plt.xlim(0,0.201)
			plt.plot(beta2, energy2, color = 'black', ls = '--', label = 'Diagonalization', linewidth=2)
	if (YLabel == "Potential"):
		plt.ylabel(r'$\mathrm{V_{0}}$ (Kelvin)', fontsize = font)
		plt.xlim(0,0.201)
		plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)

	if (YLabel == "Rotational"):
		plt.ylabel(r'$\mathrm{K_{0}^{Rot}}$ (Kelvin)', fontsize = font)
		plt.xlim(0,0.201)
		plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)

	plt.legend(loc='upper right', shadow=True, fontsize = font)
	if (variableName == "tau"):
		plt.xlim(0.000,0.00505)
		x = [0.000, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005]
		labels = ['0.000', '', '0.001', '', '0.002', ' ', '0.003', '', '0.004', ' ', '0.005']
		plt.xticks(x, labels, rotation='horizontal')


def fitFunc(var, a, b):
	return a + b*var*var

def plotfitting(val1, val2, val3):
	markersize_fig = 10
	xdata = np.linspace(0, 0.006, 10000)
	fitParams, fitCovariances = curve_fit(fitFunc, val1, val2, sigma = val3)
	plt.plot(xdata, fitFunc(xdata, fitParams[0], fitParams[1]), linestyle = '-', color = 'r', label = 'Fit', lw = 2)

def FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA):
	aa=[]
	bb=[]
	cc=[]
	dd=[]
	ee=[]
	for i in range(2,12):
		numbmolecules = i
		aa.append(numbmolecules)
		nparticle = np.array(aa)
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)
		FileToBePlot   	  = FilePlotName.SaveEnergy+".txt"
		if os.path.isfile(FileToBePlot):
			col_beads, col_tau, col_tot, err_col_tot = loadtxt(FileToBePlot,unpack=True, usecols=[0,1,5,9])
			if (col_beads.size == 1):
				bb.append(col_beads)
				cc.append(col_tau)
				dd.append(col_tot)
				ee.append(err_col_tot)
			else:
				for j in range(0,col_beads.size):
					if (col_beads[j] == 61):
						bb.append(col_beads[j])
						cc.append(col_tau[j])
						dd.append(col_tot[j])
						ee.append(err_col_tot[j])

	nbeads = np.array(bb)
	ntau   = np.array(cc)
	ntot   = np.array(dd)
	nerrtot= np.array(ee)
			
	FilePlot = FilePlotName.SaveChemPot+".pdf"
	outfile  = FilePlot

	fig = plt.figure(figsize=(6, 4), dpi=200)

	TypePlot1 = 3
	font=12
	#plt.grid(True)
	plt.xlabel('N')

	if (TypePlot1 == 1):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (Kelvin)}$', fontsize = font)
	if (TypePlot1 == 2):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle/N \ (Kelvin)}$', fontsize = font)
	if (TypePlot1 == 3):
		plt.ylabel(r'$\mathrm{\mu \ (Kelvin)}$', fontsize = font)

	plt.xlabel(r'$\mathrm{N}$', fontsize = font)

	if (TypePlot1 == 1):
		NumbRotors1  = nparticle
		TotalEnergy1 = ntot
		Error1       = nerrtot

	if (TypePlot1 == 2):
		NumbRotors1  = nparticle
		TotalEnergy1 = ntot/nparticle
		Error1       = nerrtot/nparticle

	if (TypePlot1 == 3):
		mu1 = []
		errormu1 = []
		num1 = []
    	for i in range(1,len(nbeads)):
        	ii = i-1
        	mu1.append(ntot[i] - ntot[ii])
        	num1.append(nparticle[i])
        	errormu1.append(sqrt(nerrtot[i]*nerrtot[i]+nerrtot[ii]*nerrtot[ii]))

		NumbRotors1  = num1
		TotalEnergy1 = mu1
		Error1       = errormu1

	srcfile2         = "ResultsOfPIGS/chemical_potential_unscreened.dat"
	data2            = loadtxt(srcfile2,unpack=True, usecols=[0,1])
	xdata, ydata = data2
	plt.plot(xdata, ydata, linestyle = '--', color = 'r', label = 'Dmitri', lw = 2)
	plt.errorbar(NumbRotors1, TotalEnergy1, yerr=Error1, color = 'b', ls = '-', label = 'PIGS', linewidth=2)
	plt.subplots_adjust(top=0.95, bottom=0.15, left=0.20, right=0.98, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.60, 0.98), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi = 200, format = 'pdf')
	call(["open", outfile])
	#call(["okular", outfile])
