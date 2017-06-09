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
#molecule            = "HF-C60"                                                    #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                        #change param1
molecule_rot        = "HF"
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
numbblocks	        = 5000                                                        #change param2
numbmolecules       = 4                                                            #change param3
beta     	        = 0.128                                                        #change param4
Rpt                 = 10.0                                                        #change param7
dipolemoment        = 1.86

status              = "submission"                                                 #change param9
#status              = "analysis"                                                   #change param10
status_rhomat       = "Yes"                                                        #change param10 
RUNDIR              = "scratch"

nrange              = 8 		  						                           #change param12

if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6

if (molecule_rot == "HF"):
	#step           = [0.7,1.4,2.3,4.2,7.8,5.0,2.5,1.5,0.2]  # 2 HF beta 0.512 K-1 #change param6 for 10 Angstrom
	#step           = [0.7,3.0,5.0,8.5,5.0,3.0,1.6,1.0,0.2]  # 2 HF beta 0.256 K-1 #change param6 for 10 Angstrom
	#step           = [0.7,7.0,9.5,5.5,3.0,1.5,1.0,0.7,0.2]  # 2 HF beta 0.128 K-1 #change param6 for 10 Angstrom
	step           = [0.7,1.0,1.0,1.0,1.0,0.7,0.5,0.3,0.2]  # 2 HF beta 0.032 K-1 #change param6 for 10 Angstrom
	#step            = [0.7, 15.0, 13.0, 8.0, 2.5,1.5,1.0,0.7,0.2]  # 2 HF beta 0.128 K-1 #change param6  for 10.05 Angstrom   DM 0.45
	#step            = [0.7, 10.0, 10.0, 6.0, 5.0, 2.7, 1.8, 1.2, 0.2]  # 2 HF beta 0.256 K-1 #change param6  for 10.05 Angstrom DM 0.45 
	#step            = [3.0,1.0,1.5,3.0,3.0,3.0,2.3,1.5,0.2] # 1 HF beta 0.512 K-1 #change param6 for 10 Angstrom
	#step            = [3.0,2.0,3.0,3.0,3.0,2.5,1.5,1.1,2.0] # 1 HF beta 0.256 K-1 #change param6 for 10 Angstrom
	#step            = [3.0,3.0,3.0,3.0,3.0,1.8,1.1,0.8,2.0] # 1 HF beta 0.128 K-1 Rpt = 10.05 #change param6
	#step            = [1.0, 0.03, 0.05, 0.08, 0.17, 0.25, 0.3, 0.3, 2.0] # 1 HF beta 0.128 K-1 Rpt = 2.0 #change param6
	#step            = [1.0, 0.5, 0.8, 1.0, 1.2, 1.2, 0.85, 0.6, 2.0] # 1 HF beta 0.128 K-1 Rpt = 5.0 #change param6
	#step            = [3.0,4.0,4.5,4.0,3.0,1.8,1.1,0.8,2.0] # 1 HF beta 0.128 K-1 Rpt = 10.0 #change param6

file1_name           = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
file1_name          += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path       = "/scratch/tapas/linear_rotors/" 
		final_path      = "/work/tapas/linear_rotors/"                                 #change param17

temperature         = 1.0/beta   
trunc               = 5000

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if status_rhomat == "Yes":
		support.compile_rotmat()

#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	file_output             = "Entropy-vs-"+str(var)+"-fixed-"
	file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

	call(["rm", file_output])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_entropy(status,var))

# Loop over jobs
for i in range(nrange):                                                  #change param13
	if (i>1):
 
		value        = pow(2,i) + value_min

		numbbeads    = support.dropzeros(value)
	
		tau          = beta/(value-1)

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir     = dest_path + folder_run 

		if status == "submission":
			if RUNDIR != "scratch":
				os.chdir(dest_path)
				call(["rm", "-rf", folder_run])
				call(["mkdir", folder_run])
				call(["mkdir", "-p", folder_run+"/results"])
				os.chdir(src_path)


				# copy files to running folder
				src_file      = src_path + "/IhRCOMC60.xyz"
				call(["cp", src_file, dest_dir])
				call(["cp", run_file, dest_dir])

				src_file      = src_path + "/qmc_run.input"
				call(["cp", src_file, dest_dir])

				# Write submit file for the current cycle
				os.chdir(dest_dir)

			argument1     = Rpt
			level         = support.levels(numbbeads)
			step1         = step[i]
			support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)
			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'submit_'+str(i)
			fwrite        = open(fname, 'w')
	
			if RUNDIR == "scratch":
				final_dir     = final_path + folder_run 
				input_file    = "qmc_"+argument2+str(i)+".input"
				call(["rm", "-rf", final_dir, input_file])
				print input_file
				call(["mv", "qmc.input", input_file])
				fwrite.write(support.jobstring_scratch(argument2,numbbeads,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir))
			else: 
				fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))
				

			fwrite.close()
			call(["qsub", fname, ])
			if RUNDIR != "scratch":
				os.chdir(src_path)


		if status == "analysis":

			variable = tau
			try:
				fanalyze.write(support.outputstr_entropy(numbbeads,variable,dest_dir,trunc))
			except:
				pass

if status == "analysis":
	fanalyze.close()
	call(["cat",file_output])
