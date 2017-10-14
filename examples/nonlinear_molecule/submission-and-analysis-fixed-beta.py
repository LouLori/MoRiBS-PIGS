#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = "tau"
#
TransMove           = "Yes"
RotMove             = "Yes"
#
status              = "submission"                                            
#status              = "analysis"                                            
#
NameOfServer        = "nlogn"
#NameOfServer        = "graham"
NameOfPartition     = "ntapas"
#
TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
#TypeCal             = 'ENT'
#
atom                = "H2"
atom_numb           = 1
atom_stat           = "BOSE"
atom_step_trans     = [0.2, 0.2]
atom_level          = [3, 3]
#
molecule            = "H2O"
molecule_numb       = 1
molecule_stat       = "BOLTZMANN"
molecule_step_trans = [0.2, 0.2]
molecule_level      = [4, 4]
#
molecule_rot        = "H2O"
#
numbblocks	        = 1000
numbpass            = 200
#
Rpt                 = 10.05
dipolemoment        = 0.45        #J. Chern. Phys. 73(5), 2319 (1980).

status_rhomat       = "No"                                                 
status_cagepot      = "No"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 40
loopEnd             = 102
skip                = 5

preskip             = 0
postskip            = 0

ENT_TYPE 			= "SWAPTOUNSWAP"
#ENT_TYPE 			= "SWAP"
#ENT_TYPE 			= "BROKENPATH"
#ENT_TYPE 			= "REGULARPATH"
particleA           = int(numbmolecules/2)

#extra_file_name     = "end-step-value-"
extra_file_name     = ""

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = 0.2
	parameter       = beta
	temperature     = 1.0/beta   
#==================================== MCStep ===================================# 
	if (molecule_rot == "H2"):
		step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
		step        = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
		level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

	if (molecule_rot == "HF"):
		step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3]
		step        = [1.7,1.6,1.6,1.6,1.5,1.4,1.4,1.3,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21 to 51 beta 0.1
		level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

	if (molecule_rot == "H2O"):
		step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3]
		step        = [1.7,1.6,1.6,1.6,1.5,1.4,1.4,1.3,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21 to 51 beta 0.1
		level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = 0.005
	parameter       = tau
#==================================== MCStep ===================================# 
	if (molecule_rot == "H2"):
		step_trans  = [0.3 for i in range(1000)]
		step        = [1.6 for i in range(1000)]  
		level       = [1   for i in range(1000)]

	if (molecule_rot == "HF"):
		step_trans  = [0.3 for i in range(1000)]
		step        = [2.0 for i in range(1000)]  
		level       = [1   for i in range(1000)]

#==================================Generating files for submission================#
file1_name = support.GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra_file_name)
if status   == "submission":

	if (RUNDIR == "scratch") or (NameOfServer == "graham"):
		dir_run_job = "/scratch/tapas/nonlinear-molecule/" 
	else:	
		dir_run_job     = "/work/tapas/nonlinear-molecule/"

	execution_file      = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     
	support.makeexecutionfile(src_dir,TypeCal)

if (NameOfServer == "graham"):
	dir_output      = "/scratch/tapas/nonlinear-molecule/"     
else:
	dir_output      = "/work/tapas/nonlinear-molecule/"             

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if (NameOfServer == "graham"):
		dir_run_input_pimc = "/scratch/tapas/nonlinear-molecule/"+file1_name+"PIMC"
	else:
		dir_run_input_pimc = "/work/tapas/nonlinear-molecule/"+file1_name+"PIMC"
	if (os.path.isdir(dir_run_input_pimc) == False):
		call(["rm", "-rf",  dir_run_input_pimc])
		call(["mkdir", "-p", dir_run_input_pimc])
	call(["cp", execution_file, dir_run_input_pimc])
	if status_rhomat == "Yes":
		support.compile_rotmat()
	if status_cagepot == "Yes":
		support.compile_cagepot()
		support.cagepot();
		call(["mv", "hfc60.pot", dir_run_input_pimc])

if status == "analysis":
	FileAnalysis = support.GetFileNameAnalysis(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)
	
	if (TypeCal == "ENT"):
		fanalyzeEntropy      = open(FileAnalysis.SaveEntropy, "a")
		fanalyzeEntropy.write(support.fmtAverageEntropy(status,variableName,ENT_TYPE))
	if ((TypeCal == "PIMC") or (TypeCal == "PIGS")):
		fanalyzeEnergy       = open(FileAnalysis.SaveEnergy, "a")           
		fanalyzeEnergy.write(support.fmtAverageEnergy(TypeCal,status,variableName))
		fanalyzeCorr         = open(FileAnalysis.SaveCorr, "a")           
		fanalyzeCorr.write(support.fmtAverageOrientation(status,variableName))
		fanalyzeTotalCorr    = open(FileAnalysis.SaveTotalCorr, "a")           
		fanalyzeXCorr        = open(FileAnalysis.SaveXCorr, "a")           
		fanalyzeYCorr        = open(FileAnalysis.SaveYCorr, "a")           
		fanalyzeZCorr        = open(FileAnalysis.SaveZCorr, "a")           
		fanalyzeXYCorr       = open(FileAnalysis.SaveXYCorr, "a")           

if (TypeCal == "ENT"):
	numbmolecules  *= 2
	if (variableName == "tau"):
		loopStart       = 40
	if (variableName == "beta"):
		loopStart       = 80

# Loop over jobs
if (variableName == "beta"):
	list_nb = [2,4,10,14,20,24,30,34,40]
if (variableName == "tau"):
	list_nb  = [i for i in range(loopStart, loopEnd, skip)]

iStep = 0
for i in list_nb:                                                

	if (TypeCal == 'PIMC'):

		if ((i%2) == 0):
			value    = i
		else:
			vaule    = i+1

		if (variableName == "beta"):
			beta     = tau*value
			temperature = 1.0/beta
			variable = beta
		if (variableName == "tau"):
			tau      = beta/value
			variable = tau

		numbbeads    = value
		folder_run   = file1_name+str(numbbeads)

		if status   == "submission":
			support.Submission(status, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step, step_trans, level, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, dipolemoment, status_rhomat, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep)

		if status == "analysis":

			final_dir_in_work = dir_output+folder_run
			try:
				if (TypeCal == "ENT"):
					fanalyzeEntropy.write(support.GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,ENT_TYPE))
				else:
					fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeCorr.write(GetAverageOrientation(support.numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeTotalCorr.write(support.GetAverageCorrelation("TotalCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXCorr.write(support.GetAverageCorrelation("XCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeYCorr.write(support.GetAverageCorrelation("YCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeZCorr.write(support.GetAverageCorrelation("ZCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXYCorr.write(support.GetAverageCorrelation("XYCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
			except:
				pass
	else:

		if ((i % 2) != 0):
			value    = i
		else:
			value    = i+1

		if (variableName == "beta"):
			beta     = tau*(value-1)
			temperature = 1.0/beta
			variable = beta
		if (variableName == "tau"):
			tau      = beta/(value-1)
			variable = tau

		numbbeads    = value
		folder_run   = file1_name+str(numbbeads)

		if status   == "submission":
			support.Submission(status, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step, step_trans, level, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, dipolemoment, status_rhomat, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep)

		if status == "analysis":

			final_dir_in_work = dir_output+folder_run
			try:
				if (TypeCal == "ENT"):
					fanalyzeEntropy.write(support.GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,ENT_TYPE))
				else:
					fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeCorr.write(support.GetAverageOrientation(numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeTotalCorr.write(support.GetAverageCorrelation("TotalCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXCorr.write(support.GetAverageCorrelation("XCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeYCorr.write(support.GetAverageCorrelation("YCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeZCorr.write(support.GetAverageCorrelation("ZCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXYCorr.write(support.GetAverageCorrelation("XYCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
			except:
				pass
	iStep = iStep + 1

if status == "analysis":
	if (TypeCal == "ENT"):
		fanalyzeEntropy.close()
		call(["cat",FileAnalysis.SaveEntropy])
#=========================File Checking===============================#
		SavedFile = FileAnalysis.SaveEntropy
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		try:
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		except:
			pass

	if (TypeCal == "PIGS" or TypeCal == "PIMC"):
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		fanalyzeTotalCorr.close()
		fanalyzeXCorr.close()
		fanalyzeYCorr.close()
		fanalyzeZCorr.close()
		fanalyzeXYCorr.close()
		call(["cat",FileAnalysis.SaveEnergy])
		print("")
		print("")
		#call(["cat",FileAnalysis.SaveAngDOFS])
#=========================File Checking===============================#
		try:
			SavedFile = FileAnalysis.SaveEnergy
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveTotalCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveXCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveYCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveZCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
			SavedFile = FileAnalysis.SaveXYCorr
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		except:
			pass
#=================================================================================#
#
#           for file rename
##
#file1_name1 = support.GetFileNameSubmission1(TypeCal, molecule_rot, TransMove, RotMove, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra_file_name)
##================================================================================#
#
		'''
		filetobemv = "/work/tapas/nonlinear-molecule/"+file1_name+str(numbbeads)
		filemv = "/work/tapas/nonlinear-molecule/"+file1_name1+str(numbbeads)
		print(filetobemv)
		print(filemv)
		call(["mv", filetobemv, filemv])
		'''
#

