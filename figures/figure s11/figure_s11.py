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
    path = os.getcwd()+"\\figures\\figure s11\\"
    parentfolder = os.getcwd()



    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a): 
    ssx = 0.15
    ssy = 0.2
    rect = [0.01, 0.035, 0.5, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel="Area, $A$ ($\mathrm{mm}^2$)", ylabel="Max Roughness, $k_{max}$ (mm)",pad=3)
    df = fig.load_csv(path+'\\sampling_data.csv')
    failed_list = []
    sample_list = ['AT12','AT13','AT14','AT16','AT17','AT18','AT19','AT20','AT21','AT22','AT23','AT24','AT25']
    cc_excluded = 0
    colors = ['#ec008c','#F4A261','#2A9D8F','#264653','#E9C46A']
    n1_list = []
    n2_list = []
    collapsed = 0
    for i,sample in enumerate(np.array(df)):
        # print(sample)
        folder_name = sample[0]
        area = sample[1]
        width = sample[2]
        height = sample[3]
        alpha = sample[4]
        ep = sample[5]
        k = sample [6]
        N = sample[7]

        # if alpha < 1.45 and height > 16 and width < 16:
        if folder_name in ['AT24','AT23','AT19','AT22','AT25','AT14','AT13','AT16','AT17']:
            index = np.where(folder_name==np.array(['AT24','AT23','AT19','AT22','AT25','AT14','AT13','AT16','AT17']))[0][0]
            if N ==1:
                if index<=3:
                    ax1.scatter(area, k, alpha=1,color ='gray',marker='D',zorder=3,s=60)
                else:
                    ax1.scatter(area, k, alpha=1,color =fig.colors[index-4],marker='D',zorder=3,s=60)
            else:
                if index<=3:
                    ax1.scatter(area, k, alpha=1,color ='gray',marker='s',zorder=3,s=60)
                else:
                    ax1.scatter(area, k, alpha=1,color =fig.colors[index-4],marker='s',zorder=3,s=60)
                # ax1.scatter(area, k, alpha=1,color =colors[index],marker='s',zorder=3,s=60)
            
            # Label the point with the number '1'
            if index<4:
                ax1.text(
                        area - 70, k+0.08, 'G'+str(index+1),  # Adjust the position with offsets
                        fontsize=fig.fontsize_index,            # Set font size
                        color='black',            # Set font color
                        fontweight='bold'       # Optional: Make it bold
                    )
            else:
                ax1.text(
                        area - 70, k+0.08, 'G'+str(index+1),  # Adjust the position with offsets
                        fontsize=fig.fontsize_index,            # Set font size
                        color=fig.colors[index-4],            # Set font color
                        fontweight='bold'       # Optional: Make it bold
                    )

        else:
            if N == 1:
                ax1.scatter(area, k, alpha=0.12,color ='#77CDFF',marker='D',zorder=2)
                n1_list.append([area,k])
            else:
                ax1.scatter(area, k, alpha=0.05,color ='#606060',marker='s',zorder=1)
                n2_list.append([area,k])
                collapsed = collapsed + 1
        if ep != 0.1:
            if alpha == alpha:
                failed_list.append([folder_name])
            else:
                cc_excluded = cc_excluded + 1

    


    # Figure (b)
    ssx = 0.2
    ssy = 0.4
    rect = [0.45, 0, 0.5, 1]
    fig.subplots.append(rect)
    ax2 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Tared Roughness, $\overline{k}-\overline{k}_0$ (mm)", pad=3)
    # sample_list = ['AT10','AT9','AT12','AT11','AT15']
    sample_list = ['AT25','AT14','AT13','AT16','AT17']

    # Line width
    linewidth = 2

    # EXP Data with varying colors
    for i, sample in enumerate(sample_list):
        df = fig.load_csv(os.getcwd()+'\\data\\scanning_cylinderical\\'+sample+'\\roughness_displacement.csv')
        displacement = df['Displacement']
        roughness = df['Roughness']
        ax2.plot((displacement/330)*100, roughness-roughness[0],marker='o',color=fig.colors[i])

        df = fig.load_csv(os.getcwd()+'\\data\\scanning_cylinderical\\'+sample+'\\{}_SD.csv'.format(sample))
        inc = df['Displacement']
        roughness_SD = df['Roughness_SD']
        ax2.errorbar((displacement/330)*100, roughness-roughness[0], yerr=roughness_SD, fmt='none', ecolor=fig.colors[i], capsize=0,alpha=0.5)

    # FEM Data with varying colors
    for i, sample in enumerate(sample_list):
        df = fig.load_csv(os.getcwd()+'\\data\\FEM_cylinderical_PBC\\'+sample+'\\depth.csv')
        displacement =  df.iloc[:, 0]
        roughness =  df.iloc[:, 1]
        ax2.plot((displacement)*100, roughness,'--', linewidth=linewidth, color=fig.colors[i], label='EXP')
    # Set colors
    EXP_legend = plt.Line2D([0], [0], color='gray', linestyle='-',linewidth=2, label='Exp.')
    SIM_legend = plt.Line2D([0], [0], color='gray', linestyle='--',linewidth=2, label='Sim.')
    # Set legend
    ax2.legend(handles=[EXP_legend,SIM_legend], loc='lower right')

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Saving of figure
    fig.savefig(path + 'figure_s11.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure()