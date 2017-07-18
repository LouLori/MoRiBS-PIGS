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
status              = "submission"
#status              = "analysis"

#NameOfServer        = "graham"
NameOfServer        = "nlogn"
NameOfPartition     = "tapas"

#TypeCal             = 'PIGS'
TypeCal             = 'ENT'

#molecule            = "HF-C60"                                                 
molecule            = "HF"                                                     
#molecule            = "H2"                                                   
molecule_rot        = "HF"                                                   

numbblocks	        = 400000
numbmolecules       = 4
numbpass            = 500
tau                 = 0.005

Rpt                 = 10.05
dipolemoment        = 1.86

status_rhomat       = "Yes"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 3
loopEnd             = 51
skip                = 2

preskip             = 1000
postskip            = 0
particleA           = int(numbmolecules/2)

ENT_TYPE = "SWAPTOUNSWAP"
#ENT_TYPE = "SWAP"
#ENT_TYPE = "BROKENPATH"
#ENT_TYPE = "REGULARPATH"

extra_file_name     = "-Passes"+str(numbpass)+"MemCheck"
#extra_file_name     = ""

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-e0vsbeads-"+ENT_TYPE 


file2_name          = ""                                                           #change param13
argument2           = "beads"                                                      #change param14
value_min           = 1                                                            #change param15
var                 = "beta"                                                       #change param16

src_path            = os.getcwd()
if (NameOfServer == "graham"):
    dest_path       = "/scratch/tapas/linear_rotors/"
else:
    dest_path       = "/work/tapas/linear_rotors/"
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 
			

if (NameOfServer == "graham"):
	final_path      = "/scratch/tapas/linear_rotors/"     
else:
	final_path      = "/work/tapas/linear_rotors/"             

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if (NameOfServer == "graham"):
		dest_pimc = "/scratch/tapas/linear_rotors/"+file1_name+"PIMC"
	else:
		dest_pimc = "/work/tapas/linear_rotors/"+file1_name+"PIMC"
	call(["rm", "-rf",  dest_pimc])
	call(["mkdir", "-p", dest_pimc])
	call(["cp", run_file, dest_pimc])
	if status_rhomat == "Yes":
		support.compile_rotmat()

#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	if (TypeCal == "PIGS"):
		file_output      = "Energy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		call(["rm", file_output, file_output_angularDOF])

	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+"-"+ENT_TYPE+".txt"
		call(["rm", file_output])

	if (TypeCal == "ENT"):
		fanalyze             = open(file_output, "a")
		fanalyze.write(support.fmtAverageEntropy(status,var,ENT_TYPE))
	else:
		fanalyze             = open(file_output, "a")           
		fanalyze.write(support.fmtAverageEnergy(status,var))
		fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
		fanalyze_angularDOF.write(support.fmtAverageOrientation(status,var))

if (TypeCal == "ENT"):
	numbmolecules  *= 2
	loopStart       = 49

step = [1.8 for i in range(loopStart, 1000, skip)]
# Loop over jobs
for i in range(loopStart, loopEnd, skip):                                                

	if i % 2 != 0:
		value        = i
	else:
		value        = i+value_min


	numbbeads    = value #support.dropzeros(value)
	beta         = tau*(value-1)
	temperature  = 1.0/beta

	folder_run   = file1_name+str(numbbeads)+file2_name
	dest_dir     = dest_path + folder_run 

	if status   == "submission":
		support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition)


	if status == "analysis":

		variable          = beta
		try:
			if (TypeCal == "ENT"):
				fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip,ENT_TYPE))
			else:
				fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
				fanalyze_angularDOF.write(GetAverageOrientation(support.numbbeads,variable,dest_dir,preskip,postskip))
		except:
			pass

if status == "analysis":
	fanalyze.close()
	if (TypeCal != "ENT"):
		fanalyze_angularDOF.close()
	call(["cat",file_output])
	'''
	print
	print
	call(["cat",file_output_angularDOF])
	'''
