import subprocess 					# for calling abaqus
import os 							# for file management
import sys
import time 						# for timestamps
import random
import nevergrad as ng
import numpy as np  			# for basic math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sympy as sym
from sympy import *
from sympy import simplify
from sympy import Array
from ngCallback import	ngLogger		# Custom callback function to include job ID number along with solution parameter
from concurrent import futures 			# for parallelization
import csv
import ngFunctions
from pathlib import Path
import sys


UXexp = ngFunctions.extract_Data('SS1-',[1,2,3,4,5,6,7],1, 1)

# Description:
# Unit-Cell sampling program that explores the unit-cell design space to find an
# optimal geometric design for a given optimization function.


# By David Farrell, dfarrell@g.harvard.edu
def plotAnalytical(exp,sim,fig):
	colorTab = list(mcolors.TABLEAU_COLORS)
	colorTab.extend(["black","darkgreen"])
	cc = 0
	plt.figure(fig)
	leg = []
	kk = 0
	for e,s in zip(exp,sim):
		leg.append(str(kk*15)+' deg - Exp')
		leg.append(str(kk*15)+' deg - Sim')
		plt.plot(e[0],e[1],color=colorTab[cc],linestyle="-")
		plt.plot(s[0],s[1]*12,color=colorTab[cc],linestyle="--")
		cc = cc+1
	plt.legend(leg)


def executeAbaqus(E1,E2,v12,G12,G13,G23):
	X = np.hstack([E1,E2,v12,G12,G13,G23])
	sampleNumber = hash(tuple(X))

	# Path and Folder control
	path_cwd = os.getcwd()
	fileList = ['abaqusRun','abaqusFunctions','abaqusParameters']
	sampleFolder,repeatCheck,fileName = ngFunctions.CreateWorkerDirectory(sampleNumber,E1,E2,v12,G12,G13,G23,path_db=path_cwd,fileList=fileList)
	
	# Check if sample has already been explored
	if not repeatCheck:
		FNULL = open(os.devnull, 'w')
		subprocess.call('abaqus cae noGUI='+fileName,shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()

	# Check error compared to experiments

	try:
		UXabqs = ngFunctions.extract_Data('Abaqus_Uniaxial_',[0,15,30,45,60,75,90],1,1)
		CTabqs = ngFunctions.extract_Data('Abaqus_Bending_',[1,2,3,4],1,1)
	except:
		return 5E+20

	if UXabqs != UXabqs or CTabqs != CTabqs:
		print('Sample '+str(sampleNumber)+' failed.')
		return 5E+20

	columns = [1]
	expData = UXexp
	abqsData = UXabqs

	#####################################
	##### 		Compute RMSE		#####
	#####################################
	ytSum = 0		# Sum to calculate the overall RMSE error
	ccT = 0
	kk = 0
	orientation_array= ([0,1,2,3,4,5,6])
	weights = [1]

	for exp,sim in zip(expData,abqsData):
		exp_disp = exp[0]
		exp_force = exp[1]
		sim_disp = sim[0]
		sim_force = sim[1]
		exp_offset = exp_force[-1]
		for Fe,Fs in zip(exp_force,sim_force):
			ytSum = ytSum+(Fe/exp_offset-(Fs*12)/exp_offset)**2
			ccT = ccT + 1
		plotAnalytical(exp,sim,0)
		
	returnval = np.sqrt(ytSum/(ccT+1))-0.01

	ytSum = 0		# Sum to calculate the overall RMSE error
	ccT = 0
	maxL = [35,75,35,75]
	goalY = [34.925,38.1,34.925,38.1]
	goalZ = [15.875,63.5,15.875,63.5]
	for i,CT in enumerate(CTabqs):
		dispY = CT[1]
		dispZ = CT[2]
		ytSum = ytSum+((goalY[i]-dispY[-1])/maxL[i])**2
		ytSum = ytSum+((goalZ[i]-dispZ[-1])/maxL[i])**2
		ccT = ccT+2
	#print('Cantilever RMSE: '+str(np.sqrt(ytSum/(ccT+1))))
	returnval = returnval+np.sqrt(ytSum/(ccT+1))
	os.chdir(path_cwd)
	return returnval

######################
##   	Sampling    ##
######################

if __name__=='__main__':
	#########################
	## Abaqus Python Call  ##
	#########################
	# Set up logger
	timestamp = int(time.time())
	logFile = 'log_{}.txt'.format(timestamp)
	parameterlist = ['E1','E2','v12','G12','G13','G23']
	logger = ngLogger(parameterlist,logFile,delete_folder=True)
	# Set up workers: number of 'cores' used to explore design space.
	
	numWorkers =5

	parametrization = ng.p.Instrumentation(
		E1 = ng.p.Scalar(lower=0,upper = 200),
		E2 = ng.p.Scalar(lower=0,upper = 200),
		v12 = ng.p.Scalar(lower=-15,upper = 15),
		G12 = ng.p.Scalar(lower=0,upper = 200),
		G13 = ng.p.Scalar(lower=0,upper = 200),
		G23 = ng.p.Scalar(lower=0,upper = 200))

	# Setup CMA options
	optimizer = ng.optimizers.CMA(parametrization=parametrization, budget=1000, num_workers=numWorkers)
	optimizer.parametrization.register_cheap_constraint(lambda X: ngFunctions.lamRequirement(X))
	optimizer.register_callback("tell",  logger)
	print("\n\nNevergrad initialization, beginning optimization of material parameters.\n")
	os.system('cls' if os.name == 'nt' else 'clear')	# Clear terminal
	print("Beginning optimization process, iterations will be logged in: {} file\n\nOPTIMIZATION IN PROGRESS\n".format(logFile))

	with futures.ProcessPoolExecutor(max_workers=optimizer.num_workers) as executor:
		recommendation = optimizer.minimize(executeAbaqus, executor=executor, batch_mode=True)

	logger.delete_folder = True