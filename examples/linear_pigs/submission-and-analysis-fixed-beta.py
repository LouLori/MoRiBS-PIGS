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
TransMove           = "No"
RotMove             = "Yes"
#
status              = "submission"                                            
#status              = "analysis"                                            
#
NameOfServer        = "nlogn"
#NameOfServer        = "graham"
NameOfPartition     = "ntapas"
#
#TypeCal             = "PIMC"
#TypeCal             = "PIGS"
TypeCal             = "ENT"
#
#molecule            = "HFC60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"
#
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
#
numbblocks	        = 50000
numbmolecules       = 2
numbpass            = 200
#
Rpt                 = 10.05
dipolemoment        = 1.826  #J. Chern. Phys. 73(5), 2319 (1980).
dipolemoment        = 2.0*dipolemoment

status_rhomat       = "Yes"                                                 
status_cagepot      = "No"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 40
loopEnd             = 36
skip                = 10

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
		#step       = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
		#step       = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
		step        = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
		level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

	if (molecule_rot == "HF"):
		step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3]
		#step        = [1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21,25,31,35,41,45,51 for beta 0.1
		#step        = [1.7,1.4,1.1,1.0,0.9]  # beads 21, 31, 41, 51 for beta 0.1
		#step        = [1.7,1.6,1.6,1.6,1.5,1.4,1.4,1.3,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21 to 51 beta 0.2
		step        = [2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21 to 51 beta 0.2
		#step        = [1.7,1.6,1.5,1.4,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 41,51,61,71,81,91,101 forbeta 0.2
		#step        = [1.7,1.5,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 2
		#step        = [1.0,0.9,0.8,0.7,0.7]  # beads 41, 61, 81, 101 HF
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
		dir_run_job = "/scratch/tapas/linear_rotors/" 
	else:	
		dir_run_job     = "/work/tapas/linear_rotors/"

	execution_file      = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     
	support.makeexecutionfile(src_dir,TypeCal)

if (NameOfServer == "graham"):
	dir_output      = "/scratch/tapas/linear_rotors/"     
else:
	dir_output      = "/work/tapas/linear_rotors/"             

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if (NameOfServer == "graham"):
		dir_run_input_pimc = "/scratch/tapas/linear_rotors/"+file1_name+"PIMC"
	else:
		dir_run_input_pimc = "/work/tapas/linear_rotors/"+file1_name+"PIMC"
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
		fanalyzeEnergy       = open(FileAnalysis.SaveEnergy, "a")           
		fanalyzeEnergy.write(support.fmtAverageEnergy(TypeCal,status,variableName))
		fanalyzeCorr         = open(FileAnalysis.SaveCorr, "a")           
		fanalyzeCorr.write(support.fmtAverageOrientation(status,variableName))
		fanalyzeTotalCorr    = open(FileAnalysis.SaveTotalCorr, "a")           
		fanalyzeXCorr        = open(FileAnalysis.SaveXCorr, "a")           
		fanalyzeYCorr        = open(FileAnalysis.SaveYCorr, "a")           
		fanalyzeZCorr        = open(FileAnalysis.SaveZCorr, "a")           
		fanalyzeXYCorr       = open(FileAnalysis.SaveXYCorr,"a")           

	if ((TypeCal == "PIMC") or (TypeCal == "PIGS")):
		fanalyzeEnergy       = open(FileAnalysis.SaveEnergy, "a")           
		fanalyzeEnergy.write(support.fmtAverageEnergy(TypeCal,status,variableName))
		fanalyzeCorr         = open(FileAnalysis.SaveCorr, "a")           
		fanalyzeCorr.write(support.fmtAverageOrientation(status,variableName))
		fanalyzeTotalCorr    = open(FileAnalysis.SaveTotalCorr, "a")           
		fanalyzeXCorr        = open(FileAnalysis.SaveXCorr, "a")           
		fanalyzeYCorr        = open(FileAnalysis.SaveYCorr, "a")           
		fanalyzeZCorr        = open(FileAnalysis.SaveZCorr, "a")           
		fanalyzeXYCorr       = open(FileAnalysis.SaveXYCorr,"a")           

if (TypeCal == "ENT"):
	numbmolecules  *= 2
	if (variableName == "tau"):
		loopStart       = 10
	if (variableName == "beta"):
		loopStart       = 80

# Loop over jobs
if (variableName == "beta"):
	if (TypeCal == "ENT"):
		list_nb = [4,10,14,20,24,30,34,40]
		list_nb = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
	else:
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
					fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeCorr.write(support.GetAverageOrientation(numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeTotalCorr.write(support.GetAverageCorrelation("TotalCorr", numbmolecules/2,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXCorr.write(support.GetAverageCorrelation("XCorr", numbmolecules/2,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeYCorr.write(support.GetAverageCorrelation("YCorr", numbmolecules/2,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeZCorr.write(support.GetAverageCorrelation("ZCorr", numbmolecules/2,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXYCorr.write(support.GetAverageCorrelation("XYCorr", numbmolecules/2,numbbeads,variable,final_dir_in_work,preskip,postskip))
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
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		fanalyzeTotalCorr.close()
		fanalyzeXCorr.close()
		fanalyzeYCorr.close()
		fanalyzeZCorr.close()
		fanalyzeXYCorr.close()
		call(["cat",FileAnalysis.SaveEntropy])
		print("")
		print("")
		call(["cat",FileAnalysis.SaveEnergy])
#=========================File Checking===============================#
		try:
			SavedFile = FileAnalysis.SaveEntropy
			support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
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

	if (TypeCal == "PIGS" or TypeCal == "PIMC"):
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		fanalyzeTotalCorr.close()
		fanalyzeXCorr.close()
		fanalyzeYCorr.close()
		fanalyzeZCorr.close()
		fanalyzeXYCorr.close()
		call(["cat",FileAnalysis.SaveEnergy])
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
		filetobemv = "/work/tapas/linear_rotors/"+file1_name+str(numbbeads)
		filemv = "/work/tapas/linear_rotors/"+file1_name1+str(numbbeads)
		print(filetobemv)
		print(filemv)
		call(["mv", filetobemv, filemv])
		'''
#

