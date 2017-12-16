#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import inputFile
import sys
import argparse

parser = argparse.ArgumentParser(description='It is a script file, written in Python, used to submit jobs in a queue as well as analyze output data files which are generated by the jobs submitted earlier. Note: Module support.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function), step lengths for rotational and translational motions, and levels for Bisection move (see class GetStepAndLevel) as needed.')
parser.add_argument("-d", "--DipoleMoment", type=float, help="Dipole Moment of a bipolar molecule in Debye. It is a float.")
parser.add_argument("-R", "--Rpt", type=float, help="Inter molecular spacing. It is a float.")
parser.add_argument("variable", help="Name of a variable: either beta or tau. It must be a string.", choices =["tau","beta"])
parser.add_argument("job", help="Type of a job: submission of new jobs or analyzing output files. It must be a string.", choices = ["submission", "analysis"])
parser.add_argument("cal", help="Type of calculation - it is a string: a) PIMC - Finite Temperature calculation by Path Integral Monte Carlo b) PIGS - Ground State Path Integral c) ENT - Entanglement by replica algorithm based on PIGS.", choices = ["PIMC", "PIGS", "ENT"])
parser.add_argument("--scal", help="subtype of calculations - must be defined as a string in case of ENT.", default = "SWAPTOUNSWAP", choices = ["SWAPTOUNSWAP", "BROKENPATH"])
parser.add_argument("-N", help="Number of Molecules. It must be an integer.", type = int)
parser.add_argument("-Block", help="Number of Blocks. It must be an integer", type = int)
parser.add_argument("-Pass", help="Number of Passes. It must be an integer", type = int)
parser.add_argument("--MOVECOM", action="store_true", help="allows translational motions of molecules or particles.")
parser.add_argument("--ROTMOVE", action="store_true", help="allows rotational motions of molecules or particles.")
parser.add_argument("--partition", help="allows to submit jobs in a specific cpu. It is a string.", default = "ntapas")
parser.add_argument("Molecule", help="Name of molecular system.")
parser.add_argument("Rotor", help="Name of rotor. It is needed to save rotational density matrix.")
parser.add_argument("param", type=float, help="Fixed value of beta or tau.")
parser.add_argument("--preskip", type=int, help="skips # of lines from the begining of an output file. It can be needed while analysis flag is open!", default = 0)
parser.add_argument("--postskip", type=int, help="skips # of lines from the end of an output file. It can be needed while analysis flag is open!", default = 0)
args = parser.parse_args()

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = args.variable
#
TransMove           = args.MOVECOM
RotMove             = args.ROTMOVE
#
status              = args.job
#
NameOfServer        = "nlogn"
#NameOfServer        = "graham"
NameOfPartition     = args.partition
#
TypeCal             = args.cal
#
molecule            = args.Molecule
molecule_rot        = args.Rotor
#
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
#
numbblocks	        = args.Block
numbmolecules       = args.N
numbpass            = args.Pass
#
if not args.MOVECOM:
	Rpt             = args.Rpt
dipolemoment        = args.DipoleMoment
dipolemoment        = 1.0*dipolemoment
support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
#exit()

status_rhomat       = "Yes"                                                 
status_cagepot      = "No"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

preskip             = args.preskip
postskip            = args.postskip

ENT_TYPE 			= args.scal
particleA           = int(numbmolecules/2)

#extra_file_name     = "end-step-value-"
extra_file_name     = ""

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = args.param
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = args.param
	parameter       = tau

steplevel           = inputFile.GetStepAndLevel(molecule_rot,variableName)
step_COM            = steplevel.step_trans
step_rot	        = steplevel.step
level_bisection     = steplevel.level

#==================================Generating files for submission================#
file1_name = support.GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra_file_name)
if status   == "submission":

	if (RUNDIR == "scratch") or (NameOfServer == "graham"):
		dir_run_job = "/scratch/tapas/linear_rotors/" 
	else:	
		dir_run_job     = "/work/tapas/linear_rotors/"

	execution_file      = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     
	support.makeexecutionfile(src_dir,TypeCal,ENT_TYPE)

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
	
	if (preskip >= numbblocks):
		print("")
		print("Warning!!!!!!!")
		print("============================================================================")
		print("Number of Blocks = "+str(numbblocks))
		print("Number of preskip= "+str(preskip))
		print("Error message: Number of preskip data must be less than Number of Blocks")
		print("============================================================================")
		exit()
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

list_nb = inputFile.Getbeads(TypeCal, variableName)
print(list_nb)

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
			support.Submission(status, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step_rot, step_COM, level_bisection, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, dipolemoment, status_rhomat, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep)

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
			support.Submission(status, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step_rot, step_COM, level_bisection, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, dipolemoment, status_rhomat, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep)

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
		print("")
		print("")
		call(["cat",FileAnalysis.SaveCorr])
		print("")
		print("")
		call(["cat",FileAnalysis.SaveTotalCorr])
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

