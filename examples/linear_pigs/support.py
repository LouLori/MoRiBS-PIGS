#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import math

def compile_rotmat():
	path_enter_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
	os.chdir(path_enter_linden)
	call(["make", "clean"])
	call(["make"])
	path_exit_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
	os.chdir(path_exit_linden)

def jackknife(mean,data):
	ai            = [((np.sum(data) - data[j])/(len(data) - 1.0)) for j in range(len(data))]
	deviation     = ai - mean
	devsq         = np.multiply(deviation,deviation)
	error         = sqrt(np.sum(devsq)*(len(data)-1.0)/len(data))
	return error

def levels(number):
	for j in range (10):
		jj=pow(2,j)

		if jj <= (number-1):
			level = j
		else:
			break
	return level

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)

def bconstant(molecule_rot):
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	autocminverse  = 2.1947463137e+5
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	if (molecule_rot == "HF"):
		bconst	   = 20.9561                     # in cm^-1  and it is  taken from http://webbook.nist.gov/cgi/inchi?ID=C7664393&Mask=1000#Diatomic
	return bconst

def replace(string_old, string_new, file1, file2):
	'''
	This function replaces old string by new string
	'''
	f1             = open(file1, 'r')
	f2             = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()

def beads(tau,beta):
	'''
	This function determins number of beads
	'''
	numbbeads1     =beta/tau+1
	numbbeads2     = int(round(numbbeads1,0))
	if (numbbeads2 % 2 == 0):
		numbbeads2 = numbbeads2 + 1
	return numbbeads2

def jobstring(file_name,value,thread):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	command_pimc_run = "./pimc"
	omp_thread     = str(thread)

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=%s
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, omp_thread, command_pimc_run)
	return job_string

def inputstr(numbbeads,tau,temperature):
	'''
	This function gives us the exact values of the agruments
	'''
	argu1          = "%7d"   % numbbeads
	argu2          = "%20.15f" % tau
	argu3          = "%7.5f" % temperature
	output ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+"\n"
	return output

def outputstr_energy(numbbeads,tau,dest_dir,trunc,preskip):
	'''
	This function gives us the output 
	'''
	col_block, col_pot, col_tot, col_rot = genfromtxt(dest_dir+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3], skip_header=preskip, max_rows=trunc)
	print len(col_tot)
	
	mean_pot      = np.mean(col_pot)
	mean_tot      = np.mean(col_tot)
	mean_rot      = np.mean(col_rot)

	error_pot     = jackknife(mean_pot,col_pot)
	error_tot     = jackknife(mean_tot,col_tot)
	error_rot     = jackknife(mean_rot,col_rot)
	#print i, len(col_block)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_pot, mean_tot, mean_rot, error_pot, error_tot, error_rot)
	output  += "\n"
	return output

def outputstr_angle(numbbeads,tau,dest_dir,trunc,preskip):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_theta, col_costheta1, col_theta1 = genfromtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip,max_rows=trunc)

	mean_costheta  = np.mean(col_costheta)
	mean_theta     = np.mean(col_theta)
	mean_costheta1 = np.mean(col_costheta1)
	mean_theta1    = np.mean(col_theta1)

	error_costheta = jackknife(mean_costheta,col_costheta)
	error_theta    = jackknife(mean_theta,col_theta)
	error_costheta1= jackknife(mean_costheta1,col_costheta1)
	error_theta1   = jackknife(mean_theta1,col_theta1)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_costheta, mean_theta, error_costheta, error_theta, mean_costheta1, mean_theta1, error_costheta1, error_theta1)
	output  += "\n"
	return output

def outputstr_angle1(numbbeads,tau,dest_dir,trunc,preskip):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_theta, col_costheta1, col_theta1 = genfromtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,5,6,7,8], skip_header=preskip, max_rows=trunc)

	mean_costheta  = np.mean(col_costheta)
	mean_theta     = np.mean(col_theta)
	mean_costheta1 = np.mean(col_costheta1)
	mean_theta1    = np.mean(col_theta1)

	error_costheta = jackknife(mean_costheta,col_costheta)
	error_theta    = jackknife(mean_theta,col_theta)
	error_costheta1= jackknife(mean_costheta1,col_costheta1)
	error_theta1   = jackknife(mean_theta1,col_theta1)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_costheta, mean_theta, error_costheta, error_theta, mean_costheta1, mean_theta1, error_costheta1, error_theta1)
	output  += "\n"
	return output

def fmt_energy(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. Potential', 'Avg. Total', 'Avg. rotational', 'Error of Potential', 'Error of Total', 'Error of Rotational')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)')
		output    +="\n"
		output    += '{:=<155}'.format('#')
		output    +="\n"
		return output

def fmt_angle(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. CosTheta', 'Avg. Theta', 'Error of CosTheta', 'Error of Theta', 'Avg. CosTheta1', 'Avg. Theta1', 'Error of CosTheta1', 'Error of Theta1')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)')
		output    +="\n"
		output    += '{:=<195}'.format('#')
		output    +="\n"

	return output

def fmt_angle1(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. CosTheta', 'Avg. Theta', 'Error of CosTheta', 'Error of Theta', 'Avg. CosTheta1', 'Avg. Theta1', 'Error of CosTheta1', 'Error of Theta1')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)')
		output    +="\n"
		output    += '{:=<195}'.format('#')
		output    +="\n"

	return output



def modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,distance,level,step,dipolemoment,particleA):
	'''
	This function modifies parameters in qmc_run.input
	'''
	replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("numbmolecules_input", str(numbmolecules), "qmc4.input", "qmc5.input")
	replace("distance_input", str(distance), "qmc5.input", "qmc6.input")
	replace("molecule_input", str(molecule_rot), "qmc6.input", "qmc7.input")
	replace("level_input", str(level), "qmc7.input", "qmc8.input")
	replace("dstep_input", str(step), "qmc8.input", "qmc9.input")
	replace("dipolemoment_input", str(dipolemoment), "qmc9.input", "qmc10.input")
	replace("numbpass_input", str(numbpass), "qmc10.input", "qmc11.input")
	mcskip = numbbeads*numbpass
	replace("mskip_input", str(mcskip), "qmc11.input", "qmc12.input")
	replace("numbparticle_input", str(particleA), "qmc12.input", "qmc.input")
	call(["rm", "qmc2.input"])
	call(["rm", "qmc3.input"])
	call(["rm", "qmc4.input"])
	call(["rm", "qmc5.input"])
	call(["rm", "qmc6.input"])
	call(["rm", "qmc7.input"])
	call(["rm", "qmc8.input"])
	call(["rm", "qmc9.input"])
	call(["rm", "qmc10.input"])
	call(["rm", "qmc11.input"])
	call(["rm", "qmc12.input"])


def rotmat(TypeCal,molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	if (TypeCal == 'PIMC'):
		numbbeads1		= numbbeads
	else:
		numbbeads1		= numbbeads - 1
	command_linden_run = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant(molecule))+" 150000 -1"
	print command_linden_run
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])

def jobstring_scratch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dest_pimc):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dest_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	logpath        = final_dir+"/"

	input_file     = dest_pimc+"/qmc"+file_name+str(value)+".input"
	exe_file       = dest_pimc+"/pimc"
	qmcinp         = "qmc"+file_name+str(value)+".input"

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s%s.out
#PBS -e %s%s.err
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd $PBS_O_WORKDIR
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc
mv %s /work/tapas/linear_rotors
""" % (job_name, walltime, processors, logpath, job_name, logpath, job_name, omp_thread, run_dir, output_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
	return job_string

def outputstr_entropy(numbbeads,tau,dest_dir,trunc,preskip,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	if ENT_TYPE == 'BROKENPATH':
		col_block, col_NM, col_DM = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, max_rows=trunc)
		print len(col_NM)
	
		mean_NM      = np.mean(col_NM)
		mean_DM      = np.mean(col_DM)
		mean_EN      = -log(mean_NM/mean_DM)

		error_NM     = jackknife(mean_NM,col_NM)
		error_DM     = jackknife(mean_DM,col_DM)
		error_EN     = sqrt((error_DM/mean_DM)*(error_DM/mean_DM) + (error_NM/mean_NM)*(error_NM/mean_NM))

		output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_NM, mean_DM, mean_EN, error_NM, error_DM, error_EN)
		output  += "\n"

	if ENT_TYPE == "SWAP":
		col_block, col_nm, col_dm, col_TrInv = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, max_rows=trunc)
		print len(col_block)
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_TrInv   = np.mean(col_TrInv)
		purity       = 1/mean_TrInv
		mean_EN      = -log(purity)

		error_nm     = jackknife(mean_nm,col_nm)
		error_dm     = jackknife(mean_dm,col_dm)
		error_Tr     = jackknife(mean_TrInv,col_TrInv)
		error_EN     = 0

		output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_Tr, error_EN)
		output  += "\n"

	if ENT_TYPE == "REGULARPATH":
		col_block, col_nm, col_dm, col_Tr = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, max_rows=trunc)
		print len(col_Tr)
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_Tr      = np.mean(col_Tr)
		mean_EN      = -log(mean_Tr)

		error_nm     = jackknife(mean_nm,col_nm)
		error_dm     = jackknife(mean_dm,col_dm)
		error_Tr     = jackknife(mean_Tr,col_Tr)
		error_EN     = error_Tr/mean_Tr

		output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, mean_Tr, mean_EN, error_nm, error_dm, error_Tr, error_EN)
		output  += "\n"

	return output

def fmt_entropy(status,variable,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		if ENT_TYPE == 'BROKENPATH':
			output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Entropy')
		if ENT_TYPE == 'SWAP':
			output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Purity', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		if ENT_TYPE == 'REGULARPATH':
			output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Purity', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		output    +="\n"
		output    += '{:=<205}'.format('#')
		output    +="\n"
		return output

def Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA):
	if RUNDIR != "scratch":
		os.chdir(dest_path)
		call(["rm", "-rf", folder_run])
		call(["mkdir", folder_run])
		call(["mkdir", "-p", folder_run+"/results"])
		os.chdir(src_path)

		# copy files to running folder
		call(["cp", run_file, dest_dir])

		src_file      = src_path + "/qmc_run.input"
		call(["cp", src_file, dest_dir])

		# Write submit file for the current cycle
		os.chdir(dest_dir)

	argument1     = Rpt
	level         = levels(numbbeads)
	istep         = i/skip
	step1         = step[istep]
	modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment,particleA)
	if status_rhomat == "Yes":
		rotmat(TypeCal,molecule_rot,temperature,numbbeads)

	#job submission
	fname         = 'entanglement_'+str(numbmolecules)+'particleA_'+str(particleA)+'jobsubmit_'+str(i)
	fwrite        = open(fname, 'w')

	if RUNDIR == "scratch":
		os.chdir(final_path)
		call(["rm", "-rf", folder_run])
		os.chdir(src_path)
		final_dir     = final_path + folder_run 

		#mv some files to pimc folder in /work
		input_file    = "qmc"+argument2+str(i)+".input"
		call(["mv", "qmc.input", dest_pimc+"/"+input_file])
		temperature1    = "%5.3f" % temperature
		file_rotdens    = molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
		call(["mv", file_rotdens, dest_pimc])

		if RUNIN == "CPU":
			fwrite.write(jobstring_scratch_cpu(argument2,i,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir, dest_pimc, src_path))
		else:
			fwrite.write(jobstring_scratch(argument2,i,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir, dest_pimc))
	else: 
		fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))
			

	fwrite.close()

	if (RUNIN == "CPU"):
		call(["chmod", "755", fname])
		#command_pimc_run = "./"+fname + ">"+ dest_dir+"/outpimc"+str(i)+" & "
		command_pimc_run = "./"+fname + ">outpimc"+str(i)+" & "
		print command_pimc_run
		system(command_pimc_run)
	else:
		call(["qsub", fname])
		

	if RUNDIR != "scratch":
		os.chdir(src_path)

def FileOutput(status, TypeCal, var, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc,ENT_TYPE):
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
		fanalyze.write(fmt_entropy(status,var,ENT_TYPE))
	else:
		fanalyze             = open(file_output, "a")           
		fanalyze.write(fmt_energy(status,var))
		fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
		fanalyze_angularDOF.write(fmt_angle(status,var))
		fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           


def FileClose(TypeCal):
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

def jobstring_scratch_cpu(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dest_pimc, src_path):
	'''
	This function creats jobstring for #PBS script
	'''
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dest_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"

	input_file     = dest_pimc+"/qmc"+file_name+str(value)+".input"
	exe_file       = dest_pimc+"/pimc"
	qmcinp         = "qmc"+file_name+str(value)+".input"

	job_string     = """#!/bin/bash
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd %s
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc 
mv %s /work/tapas/linear_rotors
""" % (omp_thread, run_dir, output_dir, src_path, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
	return job_string

