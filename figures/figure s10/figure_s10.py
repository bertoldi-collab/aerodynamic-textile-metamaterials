import sys
import os
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage import zoom
import pickle
from matplotlib.ticker import FormatStrFormatter

# Get the path of the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent of the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, "../"))

# Get the grandparent directory
grandparent_dir = os.path.abspath(os.path.join(parent_dir, "../"))

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)

# Add the grandparent directory to sys.path
sys.path.append(grandparent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt
from dfplt import aspect_convert
from data.uniaxial_planar.experiment_class import experiment


def figure():

    figsize = [510,200]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    # OS Cwd
    path = os.getcwd()+"\\figures\\figure s10\\"
    parentfolder = os.getcwd()



    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a): 
    figsize_a = [255,200]
    rect_a = aspect_convert(0,0,figsize_a,figsize)
    ax2 = fig.plot_pdf(rect_a,path+"PDF\\figure_s10a.pdf")


    # Figure (b)
    ssx = 0.2
    ssy = 0.4
    rect = [0.45, 0, 0.5, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Line Stress, $\sigma_L$ (MPa $\\times$ t)", pad=3)
    # sample_list = ['AT10','AT9','AT12','AT11','AT15']
    sample_list = ['AT10','AT9','AT24','AT11','AT15']

    # Line width
    linewidth = 2

    # Deformation list
    # deformation_list = [18,20,17.5,20,20]
    deformation_list = [18,20,20,20,20]

    # Set colors
    EXP_legend = plt.Line2D([0], [0], color='gray', linestyle='-',linewidth=2, label='Exp.')
    SIM_legend = plt.Line2D([0], [0], color='gray', linestyle='--',linewidth=2, label='Sim.')

    # Exp. Data
    for i, sample in enumerate(sample_list):
        AT = experiment(sample, str(int(i)),1,deformation_list[i],deformation_list[i],4,StressStrain=True,Width = 75,Height=150,Thickness=1, Verbose=False, path=parentfolder+"\\data\\uniaxial_planar\\{}\\Specimen_RawData_1.csv".format(sample))
        strain = AT.strain
        stress = AT.stress
        if sample == 'AT10':
            ax1.plot(strain*100, stress, linewidth=linewidth, color=fig.colors[i], zorder=10)
        else:
            ax1.plot(strain*100, stress, linewidth=linewidth, color=fig.colors[i])
        # Save the force-diplacement values to a csv file
        AT_df = pd.DataFrame({'Displacement': AT.disp, 'Force': AT.force})
        AT_df.to_csv(path+'{}-exp.csv'.format(sample), index=False)

    # FEM Data with varying colors
    for i, sample in enumerate(sample_list):
        df = fig.load_csv(os.getcwd()+'\\data\\FEM_planar_PBC\\unconstrained_uniaxial\\'+sample+'\\stress.csv')
        if sample =='AT11':
            strain = df.iloc[:, 0]-1
            stress = df.iloc[:, 5]
        else:
            strain = df.iloc[:, 3]-1
            stress = df.iloc[:, 8]
        ax1.plot(strain*100, stress,'--', linewidth=linewidth, color=fig.colors[i], label='FEM')

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Set legend
    ax1.legend(handles=[EXP_legend,SIM_legend], loc='upper left')

    # Saving of figure
    fig.savefig(path + 'figure_s10.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure()