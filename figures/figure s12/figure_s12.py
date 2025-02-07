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
    path = os.getcwd()+"\\figures\\figure s12\\"
    parentfolder = os.getcwd()



    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a): 
    figsize_a = [255,200]
    rect_a = aspect_convert(0,0,figsize_a,figsize)
    ax2 = fig.plot_pdf(rect_a,path+"PDF\\figure_s12a.pdf")


    # Figure (b)
    ssx = 0.2
    ssy = 0.4
    rect = [0.45, 0, 0.5, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="1st P-K Stress, $P_{22}$ (MPa)", pad=3)


    # Line width
    linewidth = 4

    # Set colors
    EXP_legend = plt.Line2D([0], [0], color='gray', linestyle='-',linewidth=2, label='Analytical')
    SIM_legend = plt.Line2D([0], [0], color='#EC008C', linestyle='--',linewidth=2, label='FEM')

    # Exp. Data
    df = fig.load_csv(os.getcwd()+'\\data\\analytical_model\\Neohookean\\strain_stress_data.csv')
    strain = df.iloc[:, 0]
    stress = df.iloc[:, 1]
    ax1.plot(strain*100, stress, linewidth=linewidth, color='gray')

    # FEM Data
    df = fig.load_csv(os.getcwd()+'\\data\\FEM_planar_PBC\\unconstrained_uniaxial\\Validation\\stress.csv')
    strain = df.iloc[:, 3]-1
    stress = df.iloc[:, 8]
    ax1.plot(strain*100, stress,'--', linewidth=linewidth, color='#EC008C', label='FEM')

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Set legend
    ax1.legend(handles=[EXP_legend,SIM_legend], loc='upper left')

    # Saving of figure
    fig.savefig(path + 'figure_s12.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure()