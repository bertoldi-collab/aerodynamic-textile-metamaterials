
import ngFunctions
import subprocess 					# for calling abaqus
import os 							# for file management
import numpy as np  			# for basic math
import sys
import shutil
import csv
sys.path.append('..')
import multiprocessing

# Description:
# Unit-Cell sampling program that explores the unit-cell design space to find an
# optimal geometric design for a giv en optimization function.

# By David Farrell, dfarrell@g.harvard.edu

def abaqusCharacterization(width,height,alpha,cyl_offset,disp_load):
	X = np.hstack([width,height,alpha])
	WorkerNumber = hash(tuple(X))	
	# Setup job name
	jobName = 'job_{}'.format(WorkerNumber)
	folder_name = "Files_"+str(WorkerNumber)

	if os.path.exists(folder_name) and os.path.isdir(folder_name):
		returnval = 1
		print('Folder already exists: '+str(WorkerNumber))
		# return returnval   
		try:
			delete = False
			with open(folder_name+'/depth.csv', 'r') as file:
				reader = csv.reader(file)
				last_row = None
				for row in reader:
					last_row = row
				if float(last_row[0]) == 0.1:
					return returnval
				else:
					print('Folder not complete, Incomplete : '+str(WorkerNumber))
					delete = True
			if delete:
				shutil.rmtree(folder_name)
		except FileNotFoundError:
			print('Folder not complete, Depth.csv not found: '+str(WorkerNumber))
			shutil.rmtree(folder_name)
	# Setup a file with all the parameters that abaqus will need to change
	ngFunctions.writeParameters(jobName,WorkerNumber,width,height,alpha,1.125,cyl_offset,disp_load)

	# create and move files and change directory to specific worker directory
	# this step is important for abaqus to avoid using similar meta-files at the same time - without this step abaqus breaks 
	DEBUG = False
	fileList = ['abaqusRun']
	owd = os.getcwd()
	WorkerFolder, owd=	ngFunctions.CreateWorkerDirectory(WorkerNumber,fileList,DEBUG)

	# run Abaqus CAE using a python script
	if DEBUG == False:
		FNULL = open(os.devnull, 'w')
		subprocess.call('abaqus cae noGUI=abaqusRun_'+str(WorkerNumber)+'.py',shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()

	return WorkerNumber

def worker(task):
    return abaqusCharacterization(*task)

def run_with_multiprocessing(numberWorkers, tasks):
    # Create a pool of workers
    with multiprocessing.Pool(processes=numberWorkers) as pool:
        # Distribute the tasks among the workers
        results = pool.map(worker, tasks)
    
    # Handle the results
    for result in results:
        print('Simulation: {} has finished successfully'.format(result))
######################
## Characterization ##
######################

if __name__=='__main__':
	#########################
	## Abaqus Optimization ##
	#########################
	# Set up workers: number of 'cores' used to explore design space.
	numWorkers = 3
	disp_load = 0.10
	cyl_offset = 0.25
	tasks = []
	for width in  np.linspace(10,22,13):
		for height in np.linspace(12,43,32):
			for per_alpha in np.linspace(0.15,0.35,9):
				alpha = np.arctan((height/2)/(width*per_alpha))
				tasks.append([width,height,alpha,cyl_offset,disp_load])
				task = tasks[-1]
	run_with_multiprocessing(numWorkers, tasks)
