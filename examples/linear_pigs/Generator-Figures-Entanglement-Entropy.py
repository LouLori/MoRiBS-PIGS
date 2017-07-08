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

TypePlot = 1

#=======================
if (TypePlot == 1):
	outfile = "PIGS-Results-Entropy-vs-beta.png"
#=======================

fig = plt.figure(figsize=(6, 4), dpi=200)

#=======================
font=18
plt.grid(True)
plt.xlim(0,0.25)
plt.xlabel(r'$\mathrm{\beta \ \  Kelvin^{-1}}$', fontsize = font)

if (TypePlot == 1):
	plt.ylabel(r'$\mathrm{S_{2}}$', fontsize = font)

x      = Specification()
#=======================

srcfile1         = "Entropy-vs-beta-fixed-tau0.005Kinv-Rpt10.05Angstrom-DipoleMoment1.86Debye-Blocks400000-System2HF-ParticleA1-preskip0-postskip0-SWAPTOUNSWAP.txt"
data1            = loadtxt(srcfile1,unpack=True, usecols=[1,5,9])
beta1, entropy1, error1 = data1

srcfile2         = "pigs_entanglement.dat"
data2            = loadtxt(srcfile2,unpack=True, usecols=[0,1])
beta2, entropy2 = data2

plt.plot(beta2, entropy2, 'g^--', markersize=5, label = 'Exact: Dmitri', lw = 1)
plt.errorbar(beta1, entropy1, yerr=error1, color = 'b', ls = '-', label = 'PIGS', linewidth=1)


#========================

plt.subplots_adjust(top=0.95, bottom=0.20, left=0.15, right=0.95, hspace=0.6, wspace=1.0)
plt.legend(bbox_to_anchor=(0.50, 0.58), loc=2, borderaxespad=0.)
plt.savefig(outfile, dpi = 200, format = 'png')

#========================

plt.show()
