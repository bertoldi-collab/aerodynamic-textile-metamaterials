
from os.path import exists
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import nevergrad as ng
import numpy as np
import time
import os 							# for file management
import time 						# for timestamps
import numpy as np  			# for basic math
import nevergrad as ng 			# for gradient-free optimization routines
import shutil 					# for file management
import glob 					# for file management
import shutil 					# for file management
import glob 					# for file management
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import csv

def ngTranslator(ep,nu,t,angle_base,L_base,H_base,disp_load):
	L_base = L_base*t
	H_base = H_base*L_base

	# Scale the points based on L_base and H_base
	ep = tuple([L_base*e for e in ep])
	nu = tuple([H_base*n for n in nu])
	return list(ep),list(nu),t,angle_base,L_base,H_base,disp_load

def rotT(T, g):
    Tprime = T
    for i in range(4):
        slices = [None] * 4
        slices[i] = slice(None)
        slices *= 2
        Tprime = g[slices].T * Tprime
    return Tprime.sum(-1).sum(-1).sum(-1).sum(-1)

def analyticalLEO(D,sig,ep,rot):
	# Rotation Matrix: Where variable 'rot' in in degrees
	cc, ss = np.cos(rot* np.pi / 180.), np.sin(rot* np.pi / 180.)
	rotMat = np.array(((cc, -ss, 0), (ss, cc, 0), (0, 0, 1)))

	Dten = np.zeros((3,3,3,3))

	D_data = {  "D1111":D[0][0],
				"D1122":D[0][1],
				"D2222":D[1][1],
				"D1133":D[0][2],
				"D2233":D[1][2],
				"D3333":D[2][2],
				"D1112":D[0][3],
				"D2212":D[1][3],
				"D3312":D[2][3],
				"D1212":D[3][3],
				"D1113":D[0][4],
				"D2213":D[1][4],
				"D3313":D[2][4],
				"D1213":D[3][4],
				"D1313":D[4][4],
				"D1123":D[0][5],
				"D2223":D[1][5],
				"D3323":D[2][5],
				"D1223":D[3][5],
				"D1323":D[4][5],
				"D2323":D[5][5]}

	for j in range(1,4):
		for i in range(1,4):
			for k in range(1,4):
				for l in range(1,4):
					if "D"+str(j)+str(i)+str(k)+str(l) in D_data:
						Dten[j-1][i-1][k-1][l-1] = D_data["D"+str(j)+str(i)+str(k)+str(l)]
					if "D"+str(i)+str(j)+str(k)+str(l) in D_data:
						Dten[j-1][i-1][k-1][l-1]  = D_data["D"+str(i)+str(j)+str(k)+str(l)]
					if "D"+str(j)+str(i)+str(l)+str(k) in D_data:
						Dten[j-1][i-1][k-1][l-1]  = D_data["D"+str(j)+str(i)+str(l)+str(k)]
					if "D"+str(k)+str(l)+str(i)+str(j) in D_data:
						Dten[j-1][i-1][k-1][l-1]  = D_data["D"+str(k)+str(l)+str(i)+str(j)]
	
	Dten_rot = rotT(Dten,rotMat)

	D1111=Dten_rot[0][0][0][0]
	D1122=Dten_rot[0][0][1][1]
	D2222=Dten_rot[1][1][1][1]
	D1133=Dten_rot[0][0][2][2]
	D2233=Dten_rot[1][1][2][2]
	D3333=Dten_rot[2][2][2][2]
	D1212=Dten_rot[0][1][0][1]
	D1313=Dten_rot[0][2][0][2]
	D2323=Dten_rot[1][2][1][2]

	D = np.array([[D1111,D1122,D1133,0,0,0],
					[D1122,D2222,D2233,0,0,0],
					[D1133,D2233,D3333,0,0,0],
					[0,0,0,D1212,0,0],
					[0,0,0,0,D1313,0],
					[0,0,0,0,0,D2323]])

	Dee = np.array([])
	Def = np.array([])
	Dfe = np.array([])
	Dff = np.array([])

	new_row = 0
	for i,s in enumerate(sig):
		for j,e in enumerate(ep):
			if s != s and e == e:
				Dee = np.concatenate((Dee,np.array([D[i][j]])))
			if s != s and e != e:
				Def = np.concatenate((Def,np.array([D[i][j]])))
			if s == 0 and e == e:
				Dfe = np.concatenate((Dfe,np.array([D[i][j]])))
			if s == 0 and e != e:
				Dff = np.concatenate((Dff,np.array([D[i][j]])))

	Dee = np.reshape(Dee,(int(np.sqrt(len(Dee))),int(np.sqrt(len(Dee)))))
	Def = np.reshape(Def,newshape=(np.size(Dee,0),int(len(Def)/np.size(Dee,0))))
	Dfe = np.reshape(Def,newshape=(int(len(Dfe)/np.size(Dee,1)),np.size(Dee,1)),order='F')
	Dff = np.reshape(Dff,newshape=(np.size(D,0)-np.size(Dee,0),np.size(D,1)-np.size(Dee,1)),order='F')

	ep_e = np.array([])
	for e in ep:
		if e == e:
			ep_e = np.concatenate((ep_e,np.array([e])))

	sig_f = np.matmul((Dee-np.matmul(Def,np.matmul(np.linalg.inv(Dff),Dfe))),ep_e)
	
	return sig_f

def orthotropicRequirement(X):#
	D = np.array([[X[1]['D1111'],X[1]['D1122'],X[1]['D1133'],X[1]['D1112'],X[1]['D1113'],X[1]['D1123']],
					[X[1]['D1122'],X[1]['D2222'],X[1]['D2233'],X[1]['D2212'],X[1]['D2213'],X[1]['D2223']],
					[X[1]['D1133'],X[1]['D2233'],X[1]['D3333'],X[1]['D3312'],X[1]['D3313'],X[1]['D3323']],
					[X[1]['D1112'],X[1]['D2212'],X[1]['D3312'],X[1]['D1212'],X[1]['D1213'],X[1]['D1223']],
					[X[1]['D1113'],X[1]['D2213'],X[1]['D3313'],X[1]['D1213'],X[1]['D1313'],X[1]['D1323']],
					[X[1]['D1123'],X[1]['D2223'],X[1]['D3323'],X[1]['D1223'],X[1]['D1323'],X[1]['D2323']]])
	
	return np.all(np.linalg.eigvals(D) > 0) and not any(np.iscomplex(np.linalg.eigvals(D)))

def lamRequirement(X):
	if X[1]['E1'] > 0 and X[1]['E2'] > 0 and X[1]['G12'] > 0 and X[1]['G13'] > 0 and X[1]['G23'] > 0 and abs(X[1]['v12']) < np.sqrt(X[1]['E1']/X[1]['E2']):
		return True
	else:
		return False

def RMSE(data1,data2):
	# Take two data sets and calculate the RMS error between them, data1 will always
	# be the primary data-set, where if data1[i] does not exist in data[2], then the result
	# will be interpolated between the two
	# Get lengths between the two datasets
	data1x = data1[0]
	data1y = data1[1]
	data2x = data2[0]
	data2y = data2[1]
	
	len_data1 = len(data1x)
	len_data2 = len(data2x)

	yTotal = 0		# 
	T = len_data1	# Counter of data points
	for p1 in range(len_data1):
		if data2x.index(p1): # If value p1 exists in second set, compare the two
			p2 = data2x.index(data1x(p1))
			yTotal = yTotal + (data2y(p2)-data1y(p1))
		else: # Otherwise interpolate between 
			temp = data2x.append(data1x(p1))
			temp.sort()
			indx = temp.index(data1x(p1))

def writeParameters(jobName,WorkerNumber,ep,nu,t,angle_base,L_base,H_base,disp_load):
	# Setup a file with all the parameters that abaqus will need to change
	wfile = open('ngParameters_'+str(WorkerNumber)+'.py', 'w')
	wfile.write('import numpy as np\r\n')
	wfile.write('WorkerNumber={}\r\n'.format(WorkerNumber))
	wfile.write('jobN = "'+jobName+'"\r\n')
	wfile.write('ep={}\r\n'.format(list(ep)))
	wfile.write('nu={}\r\n'.format(list(nu)))
	wfile.write('t={}\r\n'.format(t))
	wfile.write('angle_base={}\r\n'.format(angle_base))
	wfile.write('L_base={}\r\n'.format(L_base))
	wfile.write('H_base ={}\r\n'.format(H_base))
	wfile.write('disp_load={}\r\n'.format(disp_load))
	wfile.close()

def CreateWorkerDirectory(sampleNumber,E1,E2,v12,G12,G13,G23,path_db,fileList):
	# Name of folder
	sampleFolder='Files_'+str(sampleNumber)

	for item in fileList:
		if item == 'abaqusRun':
			shutil.copy2(item+'.py',item+'_'+str(sampleNumber)+'.py')
			fileName = item+'_'+str(sampleNumber)+'.py'


	fileName = 'abaqusRun_'+str(sampleNumber)+'.py'
	# Write to python file
	with open(fileName) as f:
		lines = f.readlines()

	with open(fileName, "w") as f:
		lines[22] = 'WorkerNumber={}\r'.format(sampleNumber)
		E1,E2,v12,G12
		lines[23] = 'E1={}\r'.format(E1)
		lines[24] = 'E2={}\r'.format(E2)
		lines[25] = 'v12={}\r'.format(v12)
		lines[26] = 'G12={}\r'.format(G12)
		lines[27] = 'G13={}\r'.format(G13)
		lines[28] = 'G23={}\r'.format(G23)
		for line in lines:
			f.write(line)

	# Obtain current path
	path_cwd = os.getcwd()
	os.chdir(path_db)
	
	# Check if current folder already exists
	repeatCheck = False
	if os.path.exists(sampleFolder):
		print('Sample '+str(sampleNumber)+' has already been sampled. Ending current simulation.')
		repeatCheck = True
		os.chdir(sampleFolder)
		return sampleFolder,repeatCheck,fileName

	# Make a folder to run abaqus simulations
	os.mkdir(sampleFolder)
	# Change path to sampling folder
	os.chdir(path_cwd)
	# Copy the abaqus files to new working folder
	shutil.move(fileName,path_db+'\\'+sampleFolder)
	shutil.copy2('SS1-1.csv',sampleFolder)
	shutil.copy2('SS1-2.csv',sampleFolder)
	shutil.copy2('SS1-3.csv',sampleFolder)
	shutil.copy2('SS1-4.csv',sampleFolder)
	shutil.copy2('SS1-5.csv',sampleFolder)
	shutil.copy2('SS1-6.csv',sampleFolder)
	shutil.copy2('SS1-7.csv',sampleFolder)

	# Change dir to db
	os.chdir(path_db)
	os.chdir(sampleFolder)

	#return the new directory folder name and the old directory path
	return sampleFolder,repeatCheck,fileName

def CleanWorkerDirectory(WorkerFolder,owd,WorkerNumber):
	'''
	for file in glob.glob('job_{}.txt'.format(WorkerNumber)):
	    shutil.copy(file, owd)
	'''
	for file in glob.glob('fig_{}.png'.format(WorkerNumber)):
	    shutil.copy(file, owd)

	#os.chdir(owd)
	shutil.rmtree(WorkerFolder)

def RMSEInterpolation(trueArray,predArray,ytSum,ccT,column):
	# Calculate root-mean square error between true val and predicated val
	# Unpack variables
	cc = 0
	if column == 1:
		trueSt = trueArray[0]
		trueSs = trueArray[1]
		predSt = predArray[0]
		predSs = predArray[1]
		for tSt in trueSt:
			[lowidx,highidx] = find_nearest(predSt,tSt)										# find closest index
			if lowidx == highidx:															# If indexs same, straight calculate (i.e. avoid divide by zero)
				ytSum = ytSum + (predSs[cc]-trueSs[cc])**2//trueSs[-1]
				cc = cc+1
			else:			
				if highidx >= len(predSs):
					break
				m = (predSs[highidx]-predSs[lowidx])/(predSt[highidx]-predSt[lowidx])		# slope of interpolant
				y_intersect = predSs[highidx]-m*predSt[lowidx]								# find the y-intersect
				pSs = m*tSt+y_intersect														# linear line equation
				ytSum = ytSum + (pSs-trueSs[cc])**2/trueSs[-1]								# Sum of error squared
				cc = cc+1

	if column == 2:
		trueSt11 = trueArray[0]
		trueSt22 = trueArray[1]
		trueSs11 = trueArray[2]
		trueSs22 = trueArray[3]
		predSt11 = predArray[0]
		predSt22 = predArray[1]
		predSs11 = predArray[2]
		predSs22 = predArray[3]
		cc11 = 0
		cc22 = 0
		for (tSt11,tSt22) in zip(trueSt11,trueSt22):
			[lowidx11,highidx11] = find_nearest(predSt11,tSt11)											# find closest index
			[lowidx22,highidx22] = find_nearest(predSt22,tSt22)											# find closest index
			if highidx11 >= len(predSs11):
					break
			if lowidx11 == highidx11:																	# If indexs same, straight calculate (i.e. avoid divide by zero)
				ytSum = ytSum + (predSs11[cc11]-trueSs11[cc11])**2/trueSs11[-1]	
				cc11 = cc11+1
			elif (predSt11[highidx11] != 0) and (predSt11[lowidx11] != 0):		
				m = (predSs11[highidx11]-predSs11[lowidx11])/(predSt11[highidx11]-predSt11[lowidx11])	# slope of interpolant
				y_intersect = predSs11[highidx11]-m*predSt11[lowidx11]									# find the y-intersect
				pSs = m*tSt11+y_intersect																# linear line equation
				ytSum = ytSum + (pSs-trueSs11[cc11])**2/trueSs11[-1]									# Sum of error squared
				cc11 = cc11+1
			if highidx22 >= len(predSs22):
					break
			if lowidx22 == highidx22:																	# If indexs same, straight calculate (i.e. avoid divide by zero)
				ytSum = ytSum + (predSs22[cc22]-trueSs22[cc22])**2//trueSs22[-1]
				cc22 = cc22+1
			elif (predSt22[highidx22] != 0) and (predSt22[lowidx22] != 0):				
				m = (predSs22[highidx22]-predSs22[lowidx22])/(predSt22[highidx22]-predSt22[lowidx22])	# slope of interpolant
				y_intersect = predSs22[highidx22]-m*predSt22[lowidx22]									# find the y-intersect
				pSs = m*tSt22+y_intersect																# linear line equation
				ytSum = ytSum + (pSs-trueSs22[cc22])**2/trueSs22[-1]									# Sum of error squared
				cc22 = cc22+1
			cc = cc11+cc22

	ccT = cc
	return ytSum,ccT

def find_nearest(array, value):
	# Basic function to find the nearest indx for a given value, if there is a true
	# value in the array then both high and low indexs will be the same.
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	if array[idx] > value:
		highidx = idx
		lowidx = idx-1
	elif array[idx] < value:
		highidx = idx+1
		lowidx = idx
	else:
		highidx = idx
		lowidx = idx

	return [lowidx,highidx]

def extract_Data(jobName,orientations,L0,Area):
	# Extract the abaqus data
	expUC = []
	#simUC = []
	cc = 0
	for angle in orientations:
		dataList = []
		csvName = jobName+str(abs(int(angle)))
		with open(csvName+'.csv') as csvfile:
			readFile = csv.reader(csvfile)
			for row in readFile:
				dataList.append(row)
		dataList = np.asarray(dataList,dtype=float)
		for row in dataList:
			if len(row)==2:
				row[0] = row[0]/L0
				row[1] = row[1]/Area
			if len(row)==4:
				row[0] = row[0]/L0
				row[1] = row[1]/L0
				row[2] = row[2]/Area
				row[3] = row[3]/Area

		expUC.append(dataList.T)
		#simUC.append(run_Analytical(dataList.T,matPar,nFibers,angle))
		#expUC[cc][4] = simUC[cc][4]
		cc = cc+1
	
	return expUC

def checkIntersection(ep,nu,L_base,t,H_base):
	# Check to see if any ligament intersect
	# Generate all relavent lines

	L_base = L_base*t
	H_base = H_base*L_base

	# Scale the points based on L_base and H_base
	ep = tuple([L_base*e for e in ep])
	nu = tuple([H_base*n for n in nu])
	intersect = False
	linePoints = ([ep[1],nu[1],0,0],
					   [ep[1],nu[1],ep[2],nu[2]],
					   [ep[1],nu[1],ep[2],nu[2]+nu[0]],
					   [ep[2],nu[2],ep[1],nu[1]-nu[0]],
					   [ep[2],nu[2],ep[3],nu[3]],
					   [ep[3],nu[3],L_base,0],
					   [ep[3],nu[3],L_base,nu[0]])
	lineEquations = []
	# Find equation for each line
	for points in linePoints:
		lineEquations.append(equationLine([points[0],points[1]],[points[2],points[3]]))
	cc = 0
	for line1 in lineEquations:
		a1,c1 = -np.array(line1)
		b1 = 1
		kk = 0
		for line2 in lineEquations:
			a2,c2 = -np.array(line2)
			b2 = 1
			if line1 != line2:
				try:
					if (a1*b2-a2*b1) == 0 or (a1*b2-a2*b1) == 0:
						intersect = False
					else:
						x0,y0 = [(b1*c2-b2*c1)/(a1*b2-a2*b1),(c1*a2-c2*a1)/(a1*b2-a2*b1)]
						if isBetween([linePoints[cc][0],linePoints[cc][1]],[linePoints[cc][2],linePoints[cc][3]],[x0,y0]) and isBetween([linePoints[kk][0],linePoints[kk][1]],[linePoints[kk][2],linePoints[kk][3]],[x0,y0]):
							intersect = True
				except:
					intersect = False
			kk = kk+1
		cc = cc+1

	return intersect

def isBetween(a, b, c):
	if max(a[0],b[0]) > c[0] and min(a[0],b[0]) < c[0] and max(a[1],b[1]) > c[1] and min(a[1],b[1]) < c[1]:
		return True
	else:
		return False

def equationLine(P1,P2):
	if (P2[0]-P1[0]) != 0:
		m = (P2[1]-P1[1])/(P2[0]-P1[0])		# Evaluate slope of the line
		b = P2[1] - m*P2[0]					# Evaluate the y-intersect
	else:
		m = 0
		b = P2[1] - m*P2[0]	
	return m,b
