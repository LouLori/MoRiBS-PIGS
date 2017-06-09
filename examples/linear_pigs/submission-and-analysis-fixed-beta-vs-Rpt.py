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
molecule_rot        = "HF"                                                         #change param2
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
numbblocks	        = 20000                                                        #change param3
numbmolecules       = 2                                                            #change param4
numbbeads           = 129                                                          #change param5
beta     	        = 0.128                                                        #change param6
dipolemoment        = 0.45                                                         #change param7
dRpt                = 0.5                                                          #change param7

status              = "submission"                                                 #change param8
status              = "analysis"                                                   #change param9
status_rhomat       = "Yes"                                                        #change param10 

nrange              = 20		  						                           #change param11

temperature         = 1.0/beta   
tau                 = beta/(numbbeads-1)

if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6

if (molecule_rot == "HF"):
	#step           = [0.7,1.4,2.3,4.2,7.8,5.0,2.5,1.5,0.2]  # 2 HF beta 0.512 K-1 #change param6
	#step           = [0.7,3.0,5.0,8.5,5.0,3.0,1.6,1.0,0.2]  # 2 HF beta 0.256 K-1 #change param6
	#step           = [0.7,7.0,9.5,5.5,3.0,1.5,1.0,0.6,0.2]  # 2 HF beta 0.128 K-1 #change param6
	#step            = [3.0,1.0,1.5,3.0,3.0,3.0,2.3,1.5,0.2] # 1 HF beta 0.512 K-1 #change param6
	#step            = [3.0,2.0,3.0,3.0,3.0,2.5,1.5,1.1,2.0] # 1 HF beta 0.256 K-1 #change param6
	#step            = [3.0,3.0,3.0,3.0,3.0,1.8,1.1,0.8,2.0] # 1 HF beta 0.128 K-1 Rpt = 10.05 #change param6
	#step            = [1.0, 0.03, 0.05, 0.08, 0.17, 0.25, 0.3, 0.3, 2.0] # 1 HF beta 0.128 K-1 Rpt = 2.0 #change param6
	#step            = [1.0, 0.5, 0.8, 1.0, 1.2, 1.2, 0.85, 0.6, 2.0] # 1 HF beta 0.128 K-1 Rpt = 5.0 #change param6
	#step            = [0.1, 0.30 ,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.75, 0.8, 0.8, 0.85, 0.9, 0.9, 0.9, 0.9, 0.9] # 1 HF beta 0.128 K-1 Rpt = 10.0 #change param6
	#step            = [0.1, 0.30 ,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.75, 0.8, 0.8, 0.85, 0.9, 0.9, 0.9, 0.9, 0.9] # 2 HF beta 0.128 K-1 Rpt = 10.0 #change param6  for dipole moment 2.86
	#step            = [0.2, 0.5, 0.55 ,0.6, 0.65, 0.65, 0.65, 0.65, 0.65, 0.7, 0.75, 0.75, 0.8, 0.8, 0.85, 0.8, 0.8, 0.8, 0.8, 0.8] # 2 HF beta 0.128 K-1 Rpt = 10.0 #change param6  for DipoleMoment 0.45 Debye beta = 0.128
	step            = [0.15, 0.5, 1.00 ,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # 2 HF beta 0.128 K-1 Rpt = 10.0 #change param6  for DipoleMoment 0.45 Debye beta = 0.256

file1_name          = "beta"+str(beta)+"Kinv-tau"+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
file1_name         += "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-e0vsRpt" 

file2_name          = "Angstrom"                                                   #change param10
argument2           = "Rpt"                                                        #change param11
value_min           = 0.5                                                            #change param12
var                 = "Rpt"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

trunc               = 20000

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
	file_output          = "Energy-vs-"+str(var)+"-fixed-"
	file_output         += "beta"+str(beta)+"Kinv-tau"+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output         += "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-trunc"+str(trunc)+".txt"
	file_output_angularDOF = "AngularDOF-vs-"+str(var)+"-fixed-"
	file_output_angularDOF+= "beta"+str(beta)+"Kinv-tau"+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output_angularDOF+= "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-trunc"+str(trunc)+".txt"
	file_output_angularDOF1 = "AngularDOF-vs-"+str(var)+"-fixed-"
	file_output_angularDOF1+= "beta"+str(beta)+"Kinv-tau"+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
	file_output_angularDOF1+= "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-trunc"+str(trunc)+"-for-zdir.txt"
	call(["rm", file_output, file_output_angularDOF, file_output_angularDOF1])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_energy(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.fmt_angle(status,var))

	fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           
	fanalyze_angularDOF1.write(support.fmt_angle1(status,var))


# Loop over jobs
for i in range(nrange):                                                  #change param13
 
	value        = i*dRpt + value_min

	Rpt          = '{:2.1f}'.format(value)

	folder_run   = file1_name+str(Rpt)+file2_name
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
		support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)
		if status_rhomat == "Yes":
			support.rotmat(molecule_rot,temperature,numbbeads)
	
		#job submission
		fname         = 'submit_'+str(i)
		fwrite        = open(fname, 'w')

		fwrite.write(support.jobstring(argument2,Rpt,numbmolecules))

		fwrite.close()
		call(["qsub", fname, ])
		os.chdir(src_path)


	if status == "analysis":

		variable = value
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
