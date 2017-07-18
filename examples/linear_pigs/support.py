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
		#bconst	   = 20.561                      # in cm^-1  and it is  taken from J. Opt. Soc. Am. Vol. 57, issue 12, page 1464, year 1967
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

def GetAverageEnergy(numbbeads,tau,dest_dir,preskip,postskip):
	'''
	This function gives us the output 
	'''
	col_block, col_pot, col_tot, col_rot, col_rot1 = genfromtxt(dest_dir+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)
	print(len(col_tot))
	
	mean_pot      = np.mean(col_pot)
	mean_tot      = np.mean(col_tot)
	mean_rot      = np.mean(col_rot)
	mean_rot1     = np.mean(col_rot1)

	error_pot     = np.std(col_pot,ddof=1)/sqrt(len(col_pot))
	error_tot     = np.std(col_tot,ddof=1)/sqrt(len(col_tot))
	error_rot     = np.std(col_rot,ddof=1)/sqrt(len(col_rot))
	error_rot1    = np.std(col_rot1,ddof=1)/sqrt(len(col_rot1))

	output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, tau, mean_pot, mean_tot, mean_rot, mean_rot1, error_pot, error_tot, error_rot, error_rot1)
	output  += "\n"
	return output

def GetAverageOrientation(numbbeads,tau,dest_dir,preskip,postskip):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_compx, col_compy, col_compz = genfromtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)
	'''
	fd             = open(dest_dir+'/results/pigs_instant.dof', 'rb')
	shape          = (-1,5) 
	dataBin        = np.fromfile(file=fd, count = -1, dtype=np.float64).reshape(shape)
	col_block      = dataBin[:,0]
	col_costheta   = dataBin[:,1]
	col_compx      = dataBin[:,2]
	col_compy      = dataBin[:,3]
	col_compz      = dataBin[:,4]
	'''

	mean_costheta  = np.mean(col_costheta)
	mean_compx     = np.mean(col_compx)
	mean_compy     = np.mean(col_compy)
	mean_compz     = np.mean(col_compz)

	error_costheta = np.std(col_costheta,ddof=1)/sqrt(len(col_costheta))
	error_compx    = np.std(col_compx,ddof=1)/sqrt(len(col_compx))
	error_compy    = np.std(col_compy,ddof=1)/sqrt(len(col_compy))
	error_compz    = np.std(col_compz,ddof=1)/sqrt(len(col_compz))

	col_abscompx   = np.fabs(col_compx)
	col_abscompy   = np.fabs(col_compy)
	col_abscompz   = np.fabs(col_compz)
	mean_abscompx  = np.mean(col_abscompx)
	mean_abscompy  = np.mean(col_abscompy)
	mean_abscompz  = np.mean(col_abscompz)
	error_abscompx = np.std(col_abscompx,ddof=1)/sqrt(len(col_abscompx))
	error_abscompy = np.std(col_abscompy,ddof=1)/sqrt(len(col_abscompy))
	error_abscompz = np.std(col_abscompz,ddof=1)/sqrt(len(col_abscompz))

	output  = '{0:10d}{1:15.5f}{2:15.5f}{3:15.5f}{4:15.5f}{5:15.5f}{6:15.5f}{7:15.5f}{8:15.5f}{9:15.5f}{10:15.5f}{11:15.5f}{12:15.5f}{13:15.5f}{14:15.5f}{15:15.5f}'.format(numbbeads, tau, mean_costheta, mean_compx, mean_compy, mean_compz, mean_abscompx, mean_abscompy, mean_abscompz, error_costheta, error_compx, error_compy, error_compz, error_abscompx, error_abscompy, error_abscompz)
	output  += "\n"
	return output

def fmtAverageEnergy(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable, 'Avg. Potential', 'Avg. Total', 'Avg. rotational', 'Avg. (E - V)', 'Error of Potential', 'Error of Total', 'Error of Rotational', 'Error of (E - V)')
		output    +="\n"
		output    +="#"
		output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('', (str(unit)), 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)')
		output    +="\n"
		output    += '{0:=<200}'.format('#')
		output    +="\n"
		return output

def fmtAverageOrientation(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}{13:^15}{14:^15}{15:^15}'.format('Beads', variable, '<sum of ei.ej>', '< compx >', '< compy >', '< compz >', '< abscompx >', '< abscompy >', '< abscompz >', 'Error: ei.ej', 'Error: compx', 'Error: compy', 'Error: compz', 'Error: abscompx', 'Error: abscompy', 'Error: abscompz')
		output    +="\n"
		output    +="#"
		output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}{13:^15}{14:^15}{15:^15}'.format('', (str(unit)), '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)')
		output    +="\n"
		output    +="#"
		output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}{13:^15}{14:^15}{15:^15}'.format('(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)', '(8)', '(9)', '(10)', '(11)', '(12)', '(13)', '(14)', '(15)', '(16)')
		output    +="\n"
		output    += '{0:=<250}'.format('#')
		output    +="\n"

	return output

def GetAverageEntropy(numbbeads,tau,dest_dir,preskip,postskip,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	if ENT_TYPE == "SWAPTOUNSWAP":
		col_block, col_nm, col_dm = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
		print(len(col_block))
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		purity       = mean_nm/mean_dm
		mean_EN      = -log(purity)

		error_nm     = np.std(col_nm,ddof=1)/sqrt(len(col_block)) 
		error_dm     = np.std(col_dm,ddof=1)/sqrt(len(col_block))
		error_Tr     = abs(purity)*sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))
		error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_Tr, error_EN)
		output  += "\n"

	if ENT_TYPE == 'BROKENPATH':
		col_block, col_nm, col_dm = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
		print(len(col_nm))
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_EN      = -log(mean_nm/mean_dm)

		error_nm     = jackknife(mean_nm,col_nm)
		error_dm     = jackknife(mean_dm,col_dm)
		error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, mean_EN, error_nm, error_dm, error_EN)
		output  += "\n"

	if ENT_TYPE == "SWAP":
		col_block, col_nm, col_dm, col_TrInv = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
		print(len(col_block))
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_TrInv   = np.mean(col_TrInv)
		purity       = 1/mean_TrInv
		mean_EN      = -log(purity)

		error_nm     = jackknife(mean_nm,col_nm)
		error_dm     = jackknife(mean_dm,col_dm)
		error_Tr     = jackknife(mean_TrInv,col_TrInv)
		error_EN     = 0 #Write the proper equation

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_Tr, error_EN)
		output  += "\n"

	if ENT_TYPE == "REGULARPATH":
		col_block, col_nm, col_dm, col_Tr = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
		print(len(col_Tr))
	
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_Tr      = np.mean(col_Tr)
		mean_EN      = -log(mean_Tr)

		error_nm     = jackknife(mean_nm,col_nm)
		error_dm     = jackknife(mean_dm,col_dm)
		error_Tr     = jackknife(mean_Tr,col_Tr)
		error_EN     = error_Tr/mean_Tr

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, tau, mean_nm, mean_dm, mean_Tr, mean_EN, error_nm, error_dm, error_Tr, error_EN)
		output  += "\n"

	return output

def fmtAverageEntropy(status,variable,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		if ENT_TYPE == 'SWAPTOUNSWAP':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Purity', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		if ENT_TYPE == 'BROKENPATH':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Entropy')
		if ENT_TYPE == 'SWAP':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Purity', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		if ENT_TYPE == 'REGULARPATH':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Purity', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		output    +="\n"
		output    += '{0:=<205}'.format('#')
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
	#mcskip = numbbeads*numbpass
	mcskip = numbpass
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
	command_linden_run = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant(molecule))+" 15000 -1"
	print(command_linden_run)
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])

def jobstring_scratch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dest_pimc):
	'''
	This function creats jobstring for #PBS script
	'''
	if (thread > 16):
		thread = 16
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

def Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition):
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
	istep         = int(i/skip)
	step1         = step[istep]
	modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment,particleA)
	if status_rhomat == "Yes":
		rotmat(TypeCal,molecule_rot,temperature,numbbeads)

	#job submission
	if (TypeCal == 'PIGS'):
		fname         = 'job-pigs-'+str(i)+'-for-'+folder_run
	if (TypeCal == 'PIMC'):
		fname         = 'job-pimc-'+str(i)+'-for-'+folder_run
	if (TypeCal == 'ENT'):
		fname         = 'job-ent-'+str(i)+'-for-'+folder_run
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
			fwrite.write(jobstring_scratch_sbatch(argument2,i,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir, dest_pimc, NameOfServer))
	else: 
		fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))
			

	fwrite.close()

	if (RUNIN == "CPU"):
		call(["chmod", "755", fname])
		#command_pimc_run = "./"+fname + ">"+ dest_dir+"/outpimc"+str(i)+" & "
		command_pimc_run = "./"+fname + ">outpimc"+str(i)+" & "
		print(command_pimc_run)
		system(command_pimc_run)
	else:
		#call(["qsub", fname])
		if (NameOfPartition == 'tapas'):
			call(["sbatch", "-p", "tapas", fname])
		else:
			call(["sbatch", fname])
		

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

def jobstring_scratch_sbatch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dest_pimc, NameOfServer):
	'''
	This function creats jobstring for #SBATCH script
	'''
	if (thread > 24):
		thread = 24
	job_name       = str(file_name)+str(value)
	walltime       = "00-48:00"
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dest_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	logpath        = dest_pimc+"/"+job_name

	input_file     = dest_pimc+"/qmc"+file_name+str(value)+".input"
	exe_file       = dest_pimc+"/pimc"
	qmcinp         = "qmc"+file_name+str(value)+".input"
	if (NameOfServer == "nlogn"):
		CommandForMove = "mv "+str(run_dir)+" /work/tapas/linear_rotors"
	if (NameOfServer == "graham"):
		CommandForMove = " "

	job_string     = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=%s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1200mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
mv %s %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
./pimc
%s
""" % (job_name, logpath, walltime, omp_thread, omp_thread, run_dir, output_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, CommandForMove)
	return job_string
