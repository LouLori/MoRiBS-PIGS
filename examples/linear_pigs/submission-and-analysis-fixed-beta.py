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
numbblocks	        = 20000                                                        #change param2
numbmolecules       = 1                                                            #change param3
beta     	        = 0.128                                                        #change param4
Rpt                 = 5.0                                                        #change param7

status              = "submission"                                                 #change param9
status              = "analysis"                                                   #change param10
status_rhomat       = "Yes"                                                        #change param10 

nrange              = 8 		  						                           #change param12

if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"

if (molecule_rot == "HF"):
	#step           = [0.7,1.4,2.3,4.2,7.8,5.0,2.5,1.5,0.2]  # 2 HF beta 0.512 K-1 #change param6
	#step           = [0.7,3.0,5.0,8.5,5.0,3.0,1.6,1.0,0.2]  # 2 HF beta 0.256 K-1 #change param6
	#step           = [0.7,7.0,9.5,5.5,3.0,1.5,1.0,0.6,0.2]  # 2 HF beta 0.128 K-1 #change param6
	#step            = [3.0,1.0,1.5,3.0,3.0,3.0,2.3,1.5,0.2] # 1 HF beta 0.512 K-1 #change param6
	#step            = [3.0,2.0,3.0,3.0,3.0,2.5,1.5,1.1,2.0] # 1 HF beta 0.256 K-1 #change param6
	#step            = [3.0,3.0,3.0,3.0,3.0,1.8,1.1,0.8,2.0] # 1 HF beta 0.128 K-1 Rpt = 10.05 #change param6
	#step            = [1.0, 0.03, 0.05, 0.08, 0.17, 0.25, 0.3, 0.3, 2.0] # 1 HF beta 0.128 K-1 Rpt = 2.0 #change param6
	#step            = [1.0, 0.5, 0.8, 1.0, 1.2, 1.2, 0.85, 0.6, 2.0] # 1 HF beta 0.128 K-1 Rpt = 5.0 #change param6
	step            = [3.0,4.0,4.5,4.0,3.0,1.8,1.1,0.8,2.0] # 1 HF beta 0.128 K-1 Rpt = 10.0 #change param6
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 
	#file1_name      = "test"

file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

temperature         = 1.0/beta   
trunc               = 10000

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
	file_output          = "Energy-vs-tau-"+str(numbmolecules)+"-"+str(molecule)
	file_output         += "-fixed-beta"+str(beta)+"-blocks"+str(numbblocks)+"Rpt"+str(Rpt)+"Angstrom-trunc"+str(trunc)+".txt"  
	file_output_angularDOF = "AngularDOF-vs-tau-"+str(numbmolecules)+"-"+str(molecule)
	file_output_angularDOF+= "-fixed-beta"+str(beta)+"-blocks"+str(numbblocks)+"Rpt"+str(Rpt)+"Angstrom-trunc"+str(trunc)+".txt"
	call(["rm", file_output, file_output_angularDOF])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_energy(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.fmt_angle(status,var))



# Loop over jobs
for i in range(nrange):                                                  #change param13
	if (i>0):
 
		value        = pow(2,i) + value_min

		numbbeads    = support.dropzeros(value)
	
		tau          = beta/(value-1)

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir     = dest_path + folder_run 

		if status == "submission":
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
			support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1)
			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'submit_'+str(i)
			fwrite        = open(fname, 'w')
	
			fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))

			fwrite.close()
			call(["qsub", fname, ])
			os.chdir(src_path)


		if status == "analysis":

			variable = tau
			fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc))
			fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc))

if status == "analysis":
	fanalyze.close()
	fanalyze_angularDOF.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
