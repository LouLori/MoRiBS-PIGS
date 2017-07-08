#!/usr/bin/env python
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('--','-','-.',':'))
		self.color2               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('--','-','-.',':'))

TypePlot = 3

#=======================
if (TypePlot == 1):
	outfile = "PIGS-Results-Energy-vs-NumberOfRotors.png"
if (TypePlot == 2):
	outfile = "PIGS-Results-EnergyOverNumberOfRotors-vs-NumberOfRotors.png"
if (TypePlot == 3):
	outfile = "PIGS-Results-ChemicalPotential-vs-NumberOfRotors.png"
#=======================

fig = plt.figure(figsize=(6, 4), dpi=200)

#=======================
font=18
plt.grid(True)
plt.xlim(1.90,8.10)
plt.xlabel('N')

if (TypePlot == 1):
	plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (Kelvin)}$', fontsize = font)
if (TypePlot == 2):
	plt.ylabel(r'$\mathrm{\langle E_{0} \rangle/N \ (Kelvin)}$', fontsize = font)
if (TypePlot == 3):
	plt.ylabel(r'$\mathrm{\mu \ (Kelvin)}$', fontsize = font)

plt.xlabel(r'$\mathrm{N}$', fontsize = font)

x      = Specification()
#=======================

srcfile1         = "Fitted-data-vs-number-of-rotors-PIGS.txt"
data1            = loadtxt(srcfile1,unpack=True, usecols=[0,1,2])
NumbRotors1, TotalEnergy1, Error1 = data1

srcfile2         = "exact-diag-data-by-tom.txt"
data2            = loadtxt(srcfile2,unpack=True, usecols=[0,1])
NumbRotors2, TotalEnergy2 = data2

if (TypePlot == 2):
	TotalEnergy1 = TotalEnergy1/NumbRotors1
	Error1       = Error1/NumbRotors1
	TotalEnergy2 = TotalEnergy2/NumbRotors2

if (TypePlot == 3):
	plt.xlim(1.90,11.10)
	mu1 = []
	errormu1 = []
	num1 = []
	mu2 = []
	num2 = []
	for i in range(1,10):
		ii = i-1
		mu1.append(TotalEnergy1[i] - TotalEnergy1[ii])
		num1.append(NumbRotors1[i])
		errormu1.append(sqrt(Error1[i]*Error1[i]+Error1[ii]*Error1[ii]))
	for i in range(1,7):
		ii = i-1
		mu2.append(TotalEnergy2[i] - TotalEnergy2[ii])
		num2.append(NumbRotors2[i])
	
	NumbRotors1  = num1
	TotalEnergy1 = mu1
	Error1       = errormu1
	NumbRotors2  = num2
	TotalEnergy2 = mu2
	print errormu1

#plt.plot(NumbRotors2, TotalEnergy2, linestyle = x.line1.next(), color = x.color1.next(), marker = x.marker1.next(), label = 'PIGS', lw = 1)
plt.plot(NumbRotors2, TotalEnergy2, 'g^--', markersize=9, label = 'Exact: Tom', lw = 2)
plt.errorbar(NumbRotors1, TotalEnergy1, yerr=Error1, color = 'b', ls = '-.', label = 'PIGS-Fit', linewidth=2)


#========================

plt.subplots_adjust(top=0.95, bottom=0.15, left=0.20, right=0.98, hspace=0.6, wspace=1.0)
plt.legend(bbox_to_anchor=(0.60, 0.98), loc=2, borderaxespad=0.)
plt.savefig(outfile, dpi = 200, format = 'png')

#========================

plt.show()
