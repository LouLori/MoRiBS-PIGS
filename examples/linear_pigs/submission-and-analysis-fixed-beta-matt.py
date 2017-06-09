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
TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
#TypeCal             = 'ENT'

#molecule            = "HF-C60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"

#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()

numbblocks	        = 20000                                                      
numbmolecules       = 2
numbpass            = 10
beta     	        = 0.1                                                       

Rpt                 = 10.0                                                     
dipolemoment        = 1.86
skip                = 10

status              = "submission"                                            
status              = "analysis"                                            
status_rhomat       = "Yes"                                                 
#RUNDIR              = "work"
RUNDIR              = "scratch"

nrange              = 101 		  						                          
trunc               = numbblocks
preskip             = 2000

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "PIMC"):
	file1_name      = "PIMC-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 


if (TypeCal == "ENT"):
	numbmolecules  *= 2
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 


file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 

final_path          = "/work/tapas/linear_rotors/"                                 #change param17

temperature         = 1.0/beta   


#==================================== MCStep ===================================# 
if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6

if (molecule_rot == "HF"):
	step           = [2.0,2.0,1.5,1.0,0.5,1.5,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0]  # 2 HF beta 0.2 K-1 #change param6 for 10 Angstrom PIGS
	#step           = [0.7, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.5, 1.5, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 1.0]  # 2 HF beta 0.2 K-1 #change param6 for 10 Angstrom PIMC

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
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
		file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
	
		file_output_angularDOF1  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF1 += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"
	
		call(["rm", file_output, file_output_angularDOF, file_output_angularDOF1])
	
	
	if (TypeCal == "PIMC"):
		file_output             = "PIMC-Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF1  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF1 += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"

		call(["rm", file_output, file_output_angularDOF, file_output_angularDOF1])

	
	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
		call(["rm", file_output])
	
	
	if (TypeCal == "ENT"):
		fanalyze             = open(file_output, "a")
		fanalyze.write(support.fmt_entropy(status,var))
	else:
		fanalyze             = open(file_output, "a")           
		fanalyze.write(support.fmt_energy(status,var))
		fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
		fanalyze_angularDOF.write(support.fmt_angle(status,var))
		fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           
	#support.FileOutput(status, TypeCal, var, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc)

# Loop over jobs
list_nb = [8,16,32,64,96]
for i in list_nb:

	if (TypeCal == 'PIMC'):

		#if (i>1 and i%skip == 0):
		if (i>1):
			value = i
			tau          = beta/value
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.outputstr_entropy(numbbeads,variable,dest_dir,trunc,preskip))
					else:
						fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc,preskip))
						fanalyze_angularDOF.write(outputstr_angle(support.numbbeads,variable,dest_dir,trunc,preskip))
						fanalyze_angularDOF1.write(outputstr_angle1(support.numbbeads,variable,dest_dir,trunc,preskip))
				except:
					pass
	else:
		if (i>1 and i % skip == 0 ):

			if i % 2 != 0:
				value        = i
			else:
				value        = i+value_min
			tau          = beta/(value-1)
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.outputstr_entropy(numbbeads,variable,dest_dir,trunc,preskip))
					else:
						fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc,preskip))
						fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc,preskip))
						fanalyze_angularDOF1.write(support.outputstr_angle1(numbbeads,variable,dest_dir,trunc,preskip))
				except:
					pass

if status == "analysis":
#	support.FileClose(TypeCal)
	fanalyze.close()
	if (TypeCal != "ENT"):
		fanalyze_angularDOF.close()
		fanalyze_angularDOF1.close()
	call(["cat",file_output])
'''
	print
	print
	call(["cat",file_output_angularDOF])
	print
	print
	call(["cat",file_output_angularDOF1])
'''
