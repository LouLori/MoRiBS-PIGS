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
status              = "analysis"                                            

NameOfServer        = "nlogn"
NameOfPartition      = "ntapas"

#NameOfServer        = "graham"
#TypeCal             = 'PIMC'
TypeCal             = 'PIGS'
#TypeCal             = 'ENT'

#molecule            = "HF-C60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"

#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()

numbblocks	        = 400000
numbmolecules       = 2
numbpass            = 500
beta     	        = 0.1

Rpt                 = 10.05
dipolemoment        = 1.86

status_rhomat       = "Yes"                                                 
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 10
loopEnd             = 51
skip                = 2

preskip             = 10000
postskip            = 0

particleA           = 1
ENT_TYPE = "SWAPTOUNSWAP"
#ENT_TYPE = "SWAP"
#ENT_TYPE = "BROKENPATH"
#ENT_TYPE = "REGULARPATH"
extra_file_name     = "-Passes"+str(numbpass)
#extra_file_name     = ""

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "PIMC"):
	file1_name      = "PIMC-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"+ENT_TYPE


file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
if (NameOfServer == "graham"):
    dest_path           = "/scratch/tapas/linear_rotors/"
else:
    dest_path           = "/work/tapas/linear_rotors/"
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 
			

if (NameOfServer == "graham"):
	final_path      = "/scratch/tapas/linear_rotors/"     
else:
	final_path      = "/work/tapas/linear_rotors/"             

temperature         = 1.0/beta   

if TypeCal == "ENT":
	intvalue = 3
else:
	intvalue =1

#==================================== MCStep ===================================# 
if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6

if (molecule_rot == "HF"):
	step           	= [0.7,1.5,1.5,1.5,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,1.9,1.8,1.7,1.6,1.6,1.5,1.5,1.5,1.4,1.4,1.4,1.4,1.4,1.3,1.3,1.3,1.3,1.3]  
					# 2 HF beta 0.1 K-1 #change param6 for 10.05 Angstrom and Dipole Moment 1.86 Debye PIGS
	#step           = [0.7, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.5, 1.5, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 1.0]  
					# 2 HF beta 0.2 K-1 #change param6 for 10 Angstrom PIMC

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

if status == "analysis":
	if (TypeCal == "PIGS"):
		file_output             = "Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		
		call(["rm", file_output, file_output_angularDOF])
	
	if (TypeCal == "PIMC"):
		file_output             = "PIMC-Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		file_output_angularDOF  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		call(["rm", file_output, file_output_angularDOF])
	
	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+"-"+ENT_TYPE+".txt"

		call(["rm", file_output])
	
	
	if (TypeCal == "ENT"):
		fanalyze             = open(file_output, "a")
		fanalyze.write(support.fmtAverageEntropy(status,var,ENT_TYPE))
	else:
		fanalyze             = open(file_output, "a")           
		fanalyze.write(support.fmtAverageEnergy(status,var))
		fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
		fanalyze_angularDOF.write(support.fmtAverageOrientation(status,var))
	#support.FileOutput(status, TypeCal, var, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc)

if (TypeCal == "ENT"):
	numbmolecules  *= 2

# Loop over jobs
#list_nb = [8,16,32,64,96,128]
#skip    = 2
#for i in list_nb:

for i in range(loopStart, loopEnd, skip):                                                
	if (TypeCal == 'PIMC'):

		if (i%skip == 0):
			value = i
			tau          = beta/value
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip))
					else:
						fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
						fanalyze_angularDOF.write(GetAverageOrientation(support.numbbeads,variable,dest_dir,preskip,postskip))
				except:
					pass
	else:
		if (i % skip == 0 ):

			if ((i % 2) != 0):
				value        = i
			else:
				value        = i+value_min
			tau          = beta/(value-1)
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip,ENT_TYPE))
					else:
						fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
						fanalyze_angularDOF.write(support.GetAverageOrientation(numbbeads,variable,dest_dir,preskip,postskip))
				except:
					pass

if status == "analysis":
#	support.FileClose(TypeCal)
	fanalyze.close()
	if (TypeCal != "ENT"):
		fanalyze_angularDOF.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
