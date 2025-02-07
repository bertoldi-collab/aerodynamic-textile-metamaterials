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
from matplotlib.patches import Arc

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

def figure(point_A, ux, uy):

    figsize = [510,200]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    # OS Cwd
    path = os.getcwd()+"\\figures\\figure s18\\"
    parentfolder = os.getcwd()

    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure (a&b)
    ssx = 0
    ssy = 0
    rect = [0.1, -0.1, 0.75, 1]
    fig.subplots.append(rect)
    ax = fig.plot([rect[0]+ssx/2, rect[1]+0.15, rect[2]-ssx/2, rect[3]-ssy/2], xlabel=None, ylabel=None, pad=3)

    ax.set_xlim(-7, 7)
    ax.set_ylim(-0.1, 7)

    # Draw Cartesian axes
    ax.arrow(0, 0, 6, 0, head_width=0.2, head_length=0.2, fc='black', ec='black')  # x-axis
    ax.arrow(0, 0, -6, 0, head_width=0.2, head_length=0.2, fc='black', ec='black')  # negative x-axis
    ax.arrow(0, 0, 0, 6, head_width=0.2, head_length=0.2, fc='black', ec='black')  # y-axis

    # Labels for Cartesian axes
    ax.text(6.3, 0, 'x', fontsize=12)
    ax.text(-6.8, 0, '-x', fontsize=12)
    ax.text(0.2, 6.1, 'y', fontsize=12)

    # Define point and polar transformation
    r = np.sqrt(point_A[0]**2 + point_A[1]**2)
    Phi = np.arctan2(point_A[1], point_A[0])

    # Calculate u_r and u_Phi
    ur = (ux * np.cos(Phi) + uy * np.sin(Phi))
    uPhi = -(ux * np.sin(Phi) - uy * np.cos(Phi))  # Flipped u_Phi

    # Calculate total displacement vector
    total_ux = ux
    total_uy = uy

    # Calculate total displacement in polar coordinates
    total_ur = (total_ux * np.cos(Phi) + total_uy * np.sin(Phi))
    total_uPhi = (total_ux * np.sin(Phi) - total_uy * np.cos(Phi)) 

    # Draw polar axes
    ax.arrow(0, 0, point_A[0], point_A[1], head_width=0.15, head_length=0.15, fc='black', ec='black')  # r-direction
    ax.arrow(point_A[0], point_A[1], point_A[1]*0.7, -point_A[0]*0.7, head_width=0.1, head_length=0.1, fc='black', ec='black')  # Phi-direction

    # Labels for polar axes
    ax.text(point_A[0]-2.9, point_A[1]-1.4, "r'-direction", fontsize=12, rotation=np.degrees(Phi))
    ax.text(point_A[0]-0.9, point_A[1]+1.1, r"$\Phi$-direction", fontsize=12, rotation=np.degrees(Phi+np.pi/2))

    # Indicate Phi
    ax.text(1.2, 0.2, r"$\Phi$", fontsize=12)
    ax.plot([0, point_A[0]/3], [0, point_A[1]/3], 'k--')

    # Velocity components for point A
    ax.arrow(point_A[0], point_A[1], ux - 0.1, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')  # ux
    ax.arrow(point_A[0], point_A[1], 0, uy - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # uy

    ax.arrow(point_A[0], point_A[1], ur*np.cos(Phi) - 0.1, ur*np.sin(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # ur
    ax.arrow(point_A[0] + ur*np.cos(Phi), point_A[1] + ur*np.sin(Phi), uPhi*np.sin(Phi) - 0.1, -uPhi*np.cos(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # uPhi

    # Label velocity components for point A
    ax.text(point_A[0]+ux+0.2, point_A[1], r"$u_x^A$", fontsize=12)
    ax.text(point_A[0], point_A[1]+uy+0.2, r"$u_y^A$", fontsize=12)
    ax.text(point_A[0]+ur*np.cos(Phi)+0.2, point_A[1]+ur*np.sin(Phi), r"$u_r^A$", fontsize=12)
    ax.text(point_A[0]+uPhi*np.sin(Phi)+0.2, point_A[1]-uPhi*np.cos(Phi), r"$u_\Phi^A$", fontsize=12)

    # Label point A and add a large dot
    ax.plot(point_A[0], point_A[1], 'ko', markersize=5)
    ax.text(point_A[0]-0.12, point_A[1]-0.5, "A", fontsize=12)

    # Define point B (mirrored in the negative x-positive y quadrant)
    point_B = [-point_A[0], point_A[1]]

    # Draw polar axes for point B
    ax.arrow(0, 0, point_B[0], point_B[1], head_width=0.15, head_length=0.15, fc='black', ec='black')  # r-direction
    ax.arrow(point_B[0], point_B[1], -point_B[1]*0.7, point_B[0]*0.7, head_width=0.1, head_length=0.1, fc='black', ec='black')  # Phi-direction

    # Labels for polar axes for point B
    ax.text(point_B[0]+0.9, point_B[1]-1.5, "r'-direction", fontsize=12, rotation=np.degrees(-Phi))
    ax.text(point_B[0]+0.5, point_B[1]+1, r"$\Phi$-direction", fontsize=12, rotation=np.degrees(-Phi+np.pi/2))

    # Velocity components for point B
    ax.arrow(point_B[0], point_B[1], -ux + 0.1, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')  # ux
    ax.arrow(point_B[0], point_B[1], 0, uy - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # uy

    ax.arrow(point_B[0], point_B[1], -ur*np.cos(Phi) + 0.1, ur*np.sin(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # ur
    ax.arrow(point_B[0] - ur*np.cos(Phi), point_B[1] + ur*np.sin(Phi), -uPhi*np.sin(Phi) + 0.1, -uPhi*np.cos(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='black', ec='black')  # uPhi

    # Label velocity components for point B
    ax.text(point_B[0]-ux-1.2, point_B[1], r"$u_x^B$", fontsize=12)
    ax.text(point_B[0], point_B[1]+uy+0.2, r"$u_y^B$", fontsize=12)
    ax.text(point_B[0]-ur*np.cos(Phi)-1.2, point_B[1]+ur*np.sin(Phi), r"$u_r^B$", fontsize=12)
    ax.text(point_B[0]-uPhi*np.sin(Phi)-0.8, point_B[1]-uPhi*np.cos(Phi), r"$u_\Phi^B$", fontsize=12)

    # Label point B and add a large dot
    ax.plot(point_B[0], point_B[1], 'ko', markersize=5)
    ax.text(point_B[0]-0.1, point_B[1]-0.5, "B", fontsize=12)

    # Draw a dark grey dashed semi-circle
    arc = Arc((0, 0), 2*r, 2*r, theta1=0, theta2=180, color='darkgrey', linestyle='--', zorder=0)
    ax.add_patch(arc)

    # Plot total displacement vector in polar coordinates
    ax.arrow(point_A[0], point_A[1], total_ur*np.cos(Phi) - 0.1, total_ur*np.sin(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='red', ec='red')  # total ur
    ax.arrow(point_A[0] + total_ur*np.cos(Phi), point_A[1] + total_ur*np.sin(Phi), total_uPhi*np.sin(Phi) - 0.1, -total_uPhi*np.cos(Phi) - 0.1, head_width=0.1, head_length=0.1, fc='red', ec='red')  # total uPhi

    # Label total displacement vector components
    ax.text(point_A[0]+total_ur*np.cos(Phi)+0.2, point_A[1]+total_ur*np.sin(Phi), r"$u_r^{total}$", fontsize=12, color='red')
    ax.text(point_A[0]+total_uPhi*np.sin(Phi)+0.2, point_A[1]-total_uPhi*np.cos(Phi), r"$u_\Phi^{total}$", fontsize=12, color='red')

    # Plot resulting vector from u_x and u_y
    ax.arrow(point_A[0], point_A[1], total_ux - 0.1, total_uy - 0.1, head_width=0.1, head_length=0.1, fc='blue', ec='blue')  # resulting vector

    # Plot u_x line starting at the end-point of the u_y vector
    ax.plot([point_A[0], point_A[0] + ux], [point_A[1] + uy, point_A[1] + uy], 'k-')

    # Grid settings
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_frame_on(False)

    # Saving of figure
    fig.savefig(path + 'figure_s18.pdf', dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    point_A = [4 * np.cos(np.pi / 6), 4 * np.sin(np.pi / 6)]
    ux = 1
    uy = 2
    figure(point_A, ux, uy)