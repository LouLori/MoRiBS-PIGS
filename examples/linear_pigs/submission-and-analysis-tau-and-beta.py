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
#molecule            = "HF-C60"                                                     #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2
numbmolecules       = 2                                                            #change param4

numbblocks	        = 20000                                                        #change param3
numbbeads           = 3
dtau                = 0.0001                                                        #change param5

Rpt                 = 10.0                                                         #change param6
dipolemoment        = 1.86

status              = "submission"                                                 #change param8
status              = "analysis"                                                   #change param8
status_rhomat       = "Yes"                                                        #change param9 
RUNDIR              = "work"

nrange              = 11  			  					                           #change param10

file1_name           = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
file1_name          += "-System"+str(numbmolecules)+str(molecule)+"-e0vstau" 

file2_name          = "Kinv"                                                      #change param13
argument2           = "tau"                                                       #change param14
var                 = "tau"                                                       #change param16

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                #change param17
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                  #change param18

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 
		final_path  = "/work/tapas/linear_rotors/"                            #change param17

trunc               = 10000

#step            = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # tau 0.001 to 0.01
step            = [1.0, 0.2, 0.3, 0.4, 0.4, 0.5, 0.7, 0.8, 0.8, 0.8, 1.0, 1.0] # tau 0.0001 to 0.0009

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
	file_output             = "Energy-vs-"+str(var)+"-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
	file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
	file_output_angularDOF1  = "AngularDOF-vs-"+str(var)+"-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"
	call(["rm", file_output, file_output_angularDOF,file_output_angularDOF1])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_energy(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.fmt_angle(status,var))
	fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           
	fanalyze_angularDOF1.write(support.fmt_angle1(status,var))

# Loop over jobs
for i in range(nrange):                                                  #change param19

	if (i>0):

		tau        = i*dtau

		beta         = tau*(numbbeads - 1)
		temperature  = 1.0/beta

		folder_run   = file1_name+str(tau)+file2_name
		dest_dir     = dest_path + folder_run 

		if status   == "submission":
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
				call(["rm", "-rf", final_dir])
				input_file    = "qmc_"+argument2+str(i)+".input"
				print input_file
				call(["mv", "qmc.input", input_file])
				fwrite.write(support.jobstring_scratch(argument2, i, numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir))
			else: 
				fwrite.write(support.jobstring(argument2, i, numbmolecules))
				

			fwrite.close()
			call(["qsub", fname])
			call(["rm", fname])
			if RUNDIR != "scratch":
				os.chdir(src_path)

		if status == "analysis":

			variable          = tau
			try:
				fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc))
				fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc))
				fanalyze_angularDOF1.write(support.outputstr_angle1(numbbeads,variable,dest_dir,trunc))
			except:
				pass

if status == "analysis":
	fanalyze.close()
	fanalyze_angularDOF.close()
	fanalyze_angularDOF1.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
	print
	print
	call(["cat",file_output_angularDOF1])
