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
    path = os.getcwd()+"\\figures\\figure s14\\"
    parentfolder = os.getcwd()



    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a): 
    figsize_a = [300,200]
    rect_a = aspect_convert(0,0,figsize_a,figsize)
    ax2 = fig.plot_pdf(rect_a,path+"PDF\\figure_s14a.pdf")

    # Figure (b)
    ssx = 0.2
    ssy = 0.4
    rect = [0.55, 0, 0.4, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Line Streess, $\sigma_{L}$ (MPa \\times t)", pad=3)

    # Line width
    linewidth = 2

    sample_list = ['AT10','AT9','AT24','AT11','AT15']
    # sample_list = ['AT24']
    nu_values = []

    for i, sample in enumerate(sample_list):
        # Uniaxial 11 (x-direction)
        df = fig.load_csv(os.getcwd()+'\\data\\FEM_planar_PBC\\unconstrained_uniaxial_poisson_ratio\\{}\\poisson_ratio_sample.csv'.format(sample))
        nu_values.append(df['Poisson_Ratio_Trad'].values[-1])

    # Plot the 5 values from v21_values with the 5 colors from fig.colors
    ax1.bar(sample_list, nu_values, color=[fig.colors[i % len(fig.colors)] for i in range(len(nu_values))])
    ax1.set_xlabel('Sample')
    ax1.set_ylabel("Initial Poisson's Ratio, $\\overline{\\nu}_{0}$")

    fig.set_ticks(ax1,[],[-0.5,0,0.5])
    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Saving of figure
    fig.savefig(path + 'figure_s14.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure()