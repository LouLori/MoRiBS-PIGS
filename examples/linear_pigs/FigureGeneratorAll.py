#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import FigureGenerator
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = "tau"
#variableName        = "beta"
#
TransMove           = "No"
RotMove             = "Yes"
#
#TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
TypeCal             = 'ENT'
#
TypePlot            = "Energy"
#TypePlot            = "ChemPot"
#TypePlot            = "CorrFunc"
#
#molecule            = "HFC60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"
#
numbblocks	        = 50000
numbmolecules       = 2
numbpass            = 200
#
Rpt                 = 10.05
dipolemoment        = 1.826        #J. Chern. Phys. 73(5), 2319 (1980).

preskip             = 0
postskip            = 0

ENT_TYPE 			= "SWAPTOUNSWAP"
#ENT_TYPE 			= "SWAP"
#ENT_TYPE 			= "BROKENPATH"
#ENT_TYPE 			= "REGULARPATH"
particleA           = int(numbmolecules/2)

extra_file_name     = ""

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = 0.2
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = 0.005
	parameter       = tau

#==================================Plotting====================================#
FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)

if (TypeCal == "ENT"):
	if (variableName == "tau"):
		ExactValueFile    = "ResultsOfPIGSENT/ExactValue-Entropy-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-Dmitri.txt"
		NumberOfMolecules, ExactEntropy = np.loadtxt(ExactValueFile, usecols=(0, 1), unpack=True)
		index = numbmolecules/2 - 1
		ExactEntropyValue = ExactEntropy[index]
		FileToBePlot   	  = FilePlotName.SaveEntropy+".txt"
		FilePlot          = FilePlotName.SaveEntropy+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureENT(FileToBePlot,FilePlot,TypeCal,variableName,parameter,ExactEntropyValue,numbmolecules,molecule,Rpt,dipolemoment)
	if (variableName == "beta"):
		ExactValueFile    = "ResultsOfPIGSENT/ENT-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-Entropy-vs-beta-fixed-tau0.005Kinv-System2HF-ParticleA1-by-Dmitri.txt"
		FileToBePlot   	  = FilePlotName.SaveEntropy+".txt"
		FilePlot          = FilePlotName.SaveEntropy+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureENT(FileToBePlot,FilePlot,TypeCal,variableName,parameter,ExactValueFile,numbmolecules,molecule,Rpt,dipolemoment)

if (TypeCal == "PIGS"):
	if (TypePlot == "CorrFunc"):
		RefPoint          = 0
		TypeCorr          = "Total"
		if (TypeCorr == "Total"):
			FileToBePlot   	  = FilePlotName.SaveTotalCorr+".txt"
			FilePlot          = FilePlotName.SaveTotalCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "XCorr"
		if (TypeCorr == "XCorr"):
			FileToBePlot   	  = FilePlotName.SaveXCorr+".txt"
			FilePlot          = FilePlotName.SaveXCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "YCorr"
		if (TypeCorr == "YCorr"):
			FileToBePlot   	  = FilePlotName.SaveYCorr+".txt"
			FilePlot          = FilePlotName.SaveYCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "ZCorr"
		if (TypeCorr == "ZCorr"):
			FileToBePlot   	  = FilePlotName.SaveZCorr+".txt"
			FilePlot          = FilePlotName.SaveZCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "XYCorr"
		if (TypeCorr == "XYCorr"):
			FileToBePlot   	  = FilePlotName.SaveXYCorr+".txt"
			FilePlot          = FilePlotName.SaveXYCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)

#End plotting ---correlation

	if (TypePlot == "Energy"):
		FileToBePlot   	  = FilePlotName.SaveEnergy+".txt"
		FilePlot          = FilePlotName.SaveEnergy+".pdf"
		call(["rm", FilePlot])
		RefPoint          = 0
		FigureGenerator.FigureEnergyPIGS(FileToBePlot,FilePlot,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)

#End plotting ---energy

	if (TypePlot == "ChemPot"):
		FigureGenerator.FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)

#End plotting ---Chemical Potential
