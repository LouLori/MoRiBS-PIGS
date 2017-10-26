#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support

Units          = support.GetUnitConverter()

molecule       = "HF"
DipoleMoment   = 1.8265                    # in Debye
BConstant      = support.GetBconst(molecule)  # in wavenumber
RCOM           = 10.05                    # in Angstrom

DipoleMomentAU = DipoleMoment/Units.AuToDebye
RCOMAU         = RCOM/Units.BOHRRADIUS
BConstantAU    = BConstant/Units.AuToCmInverse
RFactor        = RCOMAU/(DipoleMomentAU*DipoleMomentAU/BConstantAU)**(1/3)
gFactor        = DipoleMomentAU*DipoleMomentAU/(RCOMAU*RCOMAU*RCOMAU)
print(gFactor)
print(RFactor)

'''
class GetUnitConverter:
	def __init__(self):
		self.BOHRRADIUS = 0.5291772108;      # angstrom
		self.HARTREE2JL = 4.359748e-18;    	# hartree to joule  conversion factor
		self.HARTREE2KL = 3.157732e+05;    	# hartree to Kelvin conversion factor
		self.CMRECIP2KL = 1.4387672;       	# cm^-1 to Kelvin conversion factor
		self.MHZ2RCM    = 3.335640952e-5;  	# MHz to cm^-1 conversion factor

		self.AuToDebye     = 1.0/0.39343;
		self.AuToCmInverse = 219474.63137;
		self.AuToKelvin    = 315777.0;
		self.KCalperMolToCmInverse = 349.75509;

		self.HBAR  = 1.05457266; 			#  (10^-34 Js)     Planck constant
		self.AMU   = 1.6605402;  			#  (10^-27 kg)     atomic mass unit
		self.K_B   = 1.380658;   			#  (10^-23 JK^-1)  Boltzmann constant
		self.WNO2K = 0.6950356; 				# conversion from CM-1 to K
'''
