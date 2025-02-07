
import os 							# for file management
import time 						# for timestamps
import shutil 					# for file management
import glob 					# for file management
import shutil 					# for file management
import glob 					# for file management


def writeParameters(jobName,WorkerNumber,width,height,alpha,t,cyl_offset,disp_load):
	# Setup a file with all the parameters that abaqus will need to change
	wfile = open('ngParameters_'+str(WorkerNumber)+'.py', 'w')
	wfile.write('import numpy as np\r\n')
	wfile.write('WorkerNumber={}\r\n'.format(WorkerNumber))
	wfile.write('jobN = "'+jobName+'"\r\n')
	wfile.write('width={}\r\n'.format(width))
	wfile.write('height={}\r\n'.format(height))
	wfile.write('alpha={}\r\n'.format(alpha))
	wfile.write('t={}\r\n'.format(t))
	wfile.write('disp_load={}\r\n'.format(disp_load))
	wfile.write('cyl_offset={}\r\n'.format(cyl_offset))

	wfile.close()

def CreateWorkerDirectory(WorkerNumber,fileList,DEBUG):
	#name of folder
	WorkerFolder='Files_'+str(WorkerNumber)

	#obtain current path
	owd = os.getcwd()
	# INSERT COMMAS HERE
	if DEBUG == False:
		#create new folder with above defined name
		if os.path.exists(WorkerFolder):
			shutil.rmtree(WorkerFolder)
		
		os.mkdir(WorkerFolder)


		for item in fileList:
			shutil.copy2(item+'.py',item+'_temp_'+str(WorkerNumber)+'.py')
			oldFile = open(item+'_temp_'+str(WorkerNumber)+'.py', 'r')
			newFile = open(item+'_'+str(WorkerNumber)+'.py', 'w+')

			# and write everything back
			for indx,line in enumerate(oldFile):
				if indx == 33 and item == 'abaqusRun':
					newFile.write('execfile(\'ngParameters_'+str(WorkerNumber)+'.py\')\n')
				else:
					newFile.write(line)
			oldFile.close()
			newFile.close()
			os.remove(item+'_temp_'+str(WorkerNumber)+'.py')
		# #move specific files for simulatiosn into new directory
		shutil.move('abaqusRun_'+str(WorkerNumber)+'.py',WorkerFolder)
		shutil.move('ngParameters_'+str(WorkerNumber)+'.py',WorkerFolder)
	
	#change current directory into new directory
	os.chdir(WorkerFolder)
	#return the new directory folder name and the old directory path
	return WorkerFolder, owd

def CleanWorkerDirectory(WorkerFolder,owd,WorkerNumber):
	'''
	for file in glob.glob('job_{}.txt'.format(WorkerNumber)):
	    shutil.copy(file, owd)
	'''
	for file in glob.glob('fig_{}.png'.format(WorkerNumber)):
	    shutil.copy(file, owd)

	#os.chdir(owd)
	try:
		shutil.rmtree(WorkerFolder)
	except:
		pass
	return

