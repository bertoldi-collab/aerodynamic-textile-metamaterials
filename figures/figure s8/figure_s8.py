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


def figure_1():

    figsize = [510,200]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    # OS Cwd
    path = os.getcwd()+"\\figures\\figure s8\\"
    parentfolder = os.getcwd()



    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure 1 (a)
    ssx = 0.2
    ssy = 0.4
    rect = [0, 0, 0.5, 1]
    fig.subplots.append(rect)
    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel="Strain, $\epsilon$ (%)", ylabel="Line Stress, $\sigma_L$ (MPa $\\times$ t)", pad=3)
    angle_list = [0,15,30,45,60,75,90]

    SS = []
    # Create a colormap from black to light gray
    cmap = ListedColormap(['black', 'dimgray', 'gray', 'darkgray', 'silver', 'lightgray'])
    norm = plt.Normalize(vmin=0, vmax=90)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Line width
    linewidth = 2

    # Set colors
    FEM_legend = plt.Line2D([0], [0], color='#EC008C', linestyle='-',linewidth=4, label='FEM')
    EXP_legend = plt.Line2D([0], [0], color='gray', linestyle='-',linewidth=4, label='Exp.')

    # Create an empty array
    stiffness_values_exp = []

    # Exp. Data
    for i, angle in enumerate(angle_list):
        SS = experiment("SS", str(int(i)),1,5,5,4,StressStrain=True,Width = 2.25*12,Height=100,Thickness=1, Verbose=False, path=parentfolder+"\\data\\uniaxial_planar\\SS-{}\\Specimen_RawData_1.csv".format(angle))
        strain = SS.strain
        stress = SS.stress
        color = cmap(norm(angle))
        ax1.plot(strain*100, stress, linewidth=linewidth, color=color)

        # Calculate stiffness using linear fit
        coeffs = np.polyfit(strain, stress, 1)
        stiffness = coeffs[0]
        stiffness_values_exp.append(stiffness)

    # Convert stiffness values to a numpy array
    stiffness_values_exp = np.array(stiffness_values_exp)

    # Save the force-diplacement values to a csv file
    SS_df = pd.DataFrame({'Displacement': SS.disp, 'Force': SS.force})
    SS_df.to_csv(path+'SS-exp.csv', index=False)

    # Create an empty array
    stiffness_values_sim = []

    # Create a colormap from #EC008C to lighter shades with a tighter span
    cmap_fem = ListedColormap(['#EC008C', '#F061A0', '#F28AB3', '#F4B3C7', '#F6DCE0'])
    norm_fem = plt.Normalize(vmin=0, vmax=90)
    sm_fem = plt.cm.ScalarMappable(cmap=cmap_fem, norm=norm_fem)
    sm_fem.set_array([])

    # FEM Data with varying colors
    for i, angle in enumerate(angle_list):
        df = fig.load_csv(os.getcwd()+'\\data\\FEM_planar\\LAMINA_woven\\Abaqus_Uniaxial_{}.csv'.format(int(angle)))
        strain = (df['Displacement']/100)
        stress = df['Force']/(2.25)
        color = cmap_fem(norm_fem(angle))
        ax1.plot(strain*100, stress, linewidth=linewidth, color=color, label='FEM')

        # Calculate stiffness using linear fit
        coeffs = np.polyfit(strain, stress, 1)
        stiffness = coeffs[0]
        stiffness_values_sim.append(stiffness)
    # Convert stiffness values to a numpy array
    stiffness_values_sim = np.array(stiffness_values_sim)

    # Add colorbar to the figure
    cbar = fig.fig.colorbar(sm, ax=ax1, orientation='vertical')
    cbar.set_label('Angle (degrees)')
    cbar.set_ticks([0, 15, 30, 45, 60, 75, 90])

    # Set ticks
    # fig.set_ticks(ax1, [0, 5], [0,0.01,0.02,0.03])

    # Set legend
    ax1.legend(handles=[FEM_legend,EXP_legend], loc='upper left')

    # Figure 1 (a)
    ssx = 0.2
    ssy = 0.45
    rect = [0.5, 0.22, 0.5, 0.8]
    fig.subplots.append([0.555, 0, 0.5, 1])
    ax_polar = fig.plot([rect[0]+ssx/2, rect[1]+0.1, rect[2]-ssx/2, rect[3]-ssy/2], pad=3,projection='polar')

    # Convert existing axes to polar projection
    ax_polar.set_theta_direction('anticlockwise')
    # ax_polar.set_theta_offset(np.pi / 2.0)
    ax_polar.set_rlabel_position(0)
    ax_polar.grid(True)


    # Number of tiles for the polar plot:
    tiles = 4

    # Define angles for the polar plot
    angles = np.linspace(0, (tiles/2)*np.pi, (len(stiffness_values_exp)) * tiles-tiles, endpoint=False)
    # Mirror the stiffness values for the polar plot
    stiffness_temp = np.concatenate([stiffness_values_exp, stiffness_values_exp[::-1][1:]])
    stiffness_polar_exp = np.concatenate([stiffness_temp, stiffness_temp[::-1][1:]])

    # Add the first data point to the end to close the circle
    angles = np.append(angles, angles[0])

    # Plot the stiffness values on the polar plot
    ax_polar.plot(angles, stiffness_polar_exp, 'o-', color='dimgray', label='Experimental Stiffness',markersize=2)

    # Define angles for the polar plot
    angles = np.linspace(0, (tiles/2)*np.pi, (len(stiffness_values_sim)) * tiles-tiles, endpoint=False)
    # Mirror the stiffness values for the polar plot
    stiffness_temp = np.concatenate([stiffness_values_sim, stiffness_values_sim[::-1][1:]])
    stiffness_polar_sim = np.concatenate([stiffness_temp, stiffness_temp[::-1][1:]])

    # Add the first data point to the end to close the circle
    angles = np.append(angles, angles[0])

    # Plot the stiffness values on the polar plot
    ax_polar.plot(angles, stiffness_polar_sim, 'o-', color='#EC008C', label='Simulation Stiffness',markersize=2)

    # Add legend and title
    ax_polar.legend(loc='lower center', bbox_to_anchor=(0.5, -0.53))

    # Add label to the polar plot
    ax_polar.set_ylabel('Line Stiffness, $K_L$ (MPa $\\times$ t)', labelpad=30)

    # Remove the outer perimeter
    ax_polar.spines['polar'].set_visible(False)

    # Set 0 degree to the horizontal axis
    # ax_polar.set_theta_offset(np.pi / 2.0)


    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Set the limit above the max value to avoid clipping
    max_stiffness = max(max(stiffness_polar_exp), max(stiffness_values_sim))
    ax_polar.set_ylim(0, max_stiffness * 1.1)

    # Saving of figure
    fig.savefig(path + 'figure_s8.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_1()