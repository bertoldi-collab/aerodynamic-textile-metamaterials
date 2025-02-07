from math import nan
import os 						# for file management
import numpy as np  			# for basic math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import ngFunctions

def loadResults(sampleNumber):
    # Name of folder
    sampleFolder='Files_'+str(sampleNumber)
    # Obtain current path
    path_cwd = os.getcwd()
    os.chdir(sampleFolder)
    UXexp = ngFunctions.extract_Data('SS1-',[1,2,3,4,5,6,7],1,1)

    UXabqs = ngFunctions.extract_Data('Abaqus_Uniaxial_',[0,15,30,45,60,75,90],1,1)

    os.chdir(path_cwd)
    expData = [UXexp]
    abqsData = [UXabqs]
    return expData,abqsData


def plotData(xData,yData,fig,legName,colorLine,line,leg,title,maxy):
	fig = plt.figure(fig)
    # ax = fig.add_subplot(111)
    # ax.plot(self.strain,self.stress)
	plt.plot(xData,yData,color=colorLine,linestyle=line)
	# ax.set_ylim([0,maxy])
	leg.extend([legName])
	plt.title(title)	
	plt.grid()
	plt.ylabel('Force (N)')
	plt.xlabel('Displacement (mm)')
	plt.legend(leg)

######################
##     Graphing     ##
######################

if __name__=='__main__':

    sampleNumber =8380194170549771215
    expData,abqsData = loadResults(sampleNumber)
    colorTab = list(mcolors.TABLEAU_COLORS)
    colorTab.extend(["black","darkgreen"])
    leg = [[],[],[],[],[]]
    title = ['Uniaxial','Biaxial','Biaxial','Biaxial','Shear']
    for i,(exp,abqs) in enumerate(zip(expData,abqsData)):
        for j,(e,a) in enumerate(zip(exp,abqs)):
            if title[i] == 'Biaxial':
                plotData(e[0],e[2],i+1,'Experimental '+str(j),colorTab[j],"-",leg[i],title[i],0)
                plotData(e[1],e[3],i+1,'Experimental '+str(j),colorTab[j],"-.",leg[i],title[i],0)
                plotData(a[0],a[2],i+1,'Simulation '+str(j),colorTab[j],"--",leg[i],title[i],0)
                plotData(a[1],a[3],i+1,'Simulation '+str(j),colorTab[j],":",leg[i],title[i],0)
            else:
                plotData(e[0],e[1],i+1,'Experimental '+str(j),colorTab[j],"-",leg[i],title[i],0)
                plotData(a[0],a[1]*12,i+1,'Simulation '+str(j),colorTab[j],"--",leg[i],title[i],0)

    plt.show()

    

