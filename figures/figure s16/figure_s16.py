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
    path = os.getcwd()+"\\figures\\figure s16\\"
    parentfolder = os.getcwd()

    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a&b)
    ssx = 0.2
    ssy = 0.4
    rect = [0.0, 0, 0.45, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Line Stress, $\sigma_L$ (MPa $\\times$ t)", pad=3)
    rect = [0.5, 0, 0.45, 1]
    fig.subplots.append(rect)
    ax2 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Line Stress, $\sigma_L$ (MPa $\\times$ t)", pad=3)

    sample_list = ['AT24']

    # Line width
    linewidth = 2

    # Deformation list
    deformation_list = [20]

    # Exp. Data
    for i, sample in enumerate(sample_list):
        AT = experiment(sample, str(int(i)),1,deformation_list[i],deformation_list[i],4,StressStrain=True,Width = 75,Height=150,Thickness=1, Verbose=False, path=parentfolder+"\\data\\uniaxial_planar\\{}\\Specimen_RawData_1.csv".format(sample))
        
        strain_unprocessed = AT.disp_unprocessed/AT.height
        stress_unprocessed = AT.force_unprocessed/(AT.width)
        
        strain = AT.strain
        stress = AT.stress
        
        ax1.plot(strain_unprocessed*100, stress_unprocessed, linewidth=linewidth, color=fig.colors[2], zorder=10)
        # Split the data into 4 segments
        num_segments = 4
        segment_length = len(strain_unprocessed) // num_segments
        colors = [fig.colors[i] for i in range(num_segments)]
        labels = ['Cycle 1', 'Cycle 2', 'Cycle 3', 'Cycle 4']

        for j in range(num_segments):
            start_idx = j * segment_length
            if j == num_segments - 1:
                end_idx = len(strain_unprocessed)
            else:
                end_idx = (j + 1) * segment_length

            ax1.plot(strain_unprocessed[start_idx:end_idx] * 100, stress_unprocessed[start_idx:end_idx], 
                 linewidth=linewidth, color=colors[j], zorder=10, label=labels[j])

        ax1.legend()
        ax2.plot(strain*100, stress, linewidth=linewidth, color='#008CEC', zorder=10, label='Cycle Average')
        ax2.legend()

        # Save the force-diplacement values to a csv file
        AT_df = pd.DataFrame({'Displacement': AT.disp, 'Force': AT.force})
        AT_df.to_csv(path+'{}-exp.csv'.format(sample), index=False)

    fig.set_ticks(ax1,[0,5,10],[0,0.5,1.0])
    fig.set_ticks(ax2,[0,5,10],[0,0.5,1.0])

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Saving of figure
    fig.savefig(path + 'figure_s16.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure()