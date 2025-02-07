import sys
import os
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter
import colorcet as cc

# Get the path of the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent of the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, "../"))

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt
from dfplt import aspect_convert

def figure_3():
    def scientific(x, pos):
        return f'{x:.0e}'
    def reynolds_to_wind_speed(x):
        return x*dynamic_viscosity/diameter_cylinder

    def wind_speed_to_reynolds(x):
        return x*diameter_cylinder/dynamic_viscosity
    scaler = 1
    figsize = [249.5,460]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    colors =['#CDCDCD','#333333']
    cm = fig.custom_cm(colors)
    # cm = fig.custom_cm(cc.CET_L1).reversed()
    cm = plt.cm.get_cmap('jet')

    # Constant:
    diameter_cylinder = 0.05715     # 2.25inch cylinder (mm)
    dynamic_viscosity = 1.516*10**-5 # m^2/s  (20 degree C)
    markersize = 2.5
    linewidth = 1.5

    # OS Cwd
    path = os.getcwd()+"\\figures\\figure 3\\"
    ax_list = []
    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    
    # figure 3 (a)
    figsize_a = [figsize[0],160]
    rect_a = aspect_convert(0,0.667,figsize_a,figsize)
    ax1 = fig.plot_pdf(rect_a,path+"figure_3_a\\PDF\\figure_3_a.pdf") 
    # Add LaTeX text to specific positions
    ax1_text = fig.plot(rect_a, xlabel=None, ylabel=None,pad=3)
    ax1_text.text(0.06, 0.20, '$U$', fontsize=fig.fontsize_label+2, color='#1B99D5')  # Adding $V$ at position (1, 10)
    ax1_text.text(0.81, 0.81, '$F$', fontsize=fig.fontsize_label+2, color='#1B99D5')    # Adding $F$ at position (3, 9)
    ax1_text.text(0.84, 0.58, r'$\epsilon$', fontsize=fig.fontsize_label+2, color='#1B99D5')  # Adding $\epsilon$ at position (5, 13)
    ax1_text.axis('off')
    ax1_text.set_xlim([0,1])
    ax1_text.set_ylim([0,1])

    df_smooth_cylinder = fig.load_csv(path+'smooth_cylinder.csv')
    Re_sc = df_smooth_cylinder['Re']
    Cd_sc = df_smooth_cylinder['Cd']

    # Figure 3 (b)
    ssx = 0.35
    ssy = 0.04
    rect = [0.05, 0.33, 1, 0.3]
    fig.subplots.append([0.0, 0.36, 1, 0.35])
    ax2 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel=None, ylabel="Drag Coefficient, $C_d$",pad=3)
    ax2.plot(Re_sc, Cd_sc, '--k', linewidth=linewidth,label='Smooth Cylinder [7]')

    sample = 'KNT'
    colors =['#CDCDCD','#333333']
    # cm = fig.custom_cm(colors)
    increments = [0,10,20,30,40,50]
    rect = patches.Rectangle((wind_speed_to_reynolds(6.5),0.6), wind_speed_to_reynolds(20-6.5), 1.4-0.6, linewidth=1, edgecolor='none', facecolor='black', alpha=0.3)
    ax2.add_patch(rect)

    for i,inc in enumerate(increments):
        df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel\\'+sample+'\\'+sample+'_'+str(inc)+'_AVG.csv')
        Re = df['Re']
        Cd = df['Cd']
        ax2.plot(Re,Cd,marker='o',color=cm(inc/increments[-1]),markersize=markersize,linewidth=linewidth)
    # Define the normalization for the color mapper
    norm = Normalize(vmin=0, vmax=increments[-1]*1e-1/0.33)
    cbar1 = fig.add_cm_plot(ax2,cm,norm,ticks=np.arange(0,15,2),label='Strain, $\\epsilon$ (%)')
    ax2.legend(fontsize=7)
    ax_list.append(ax2)

    # Figure 3 (c)
    rect = [0.05, 0.05, 1, 0.3]
    fig.subplots.append([0.00, 0.01, 1, 0.35])
    ax3 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel=None, ylabel="Drag Coefficient, $C_d$",pad=1.5)
    ax3.plot(Re_sc, Cd_sc, '--k', linewidth=linewidth,label='Smooth Cylinder [7]')
    sample = 'AT24'
    colors =['#fbcce8','#EC008C']
    # cm = fig.custom_cm(colors)
    increments = [0,3,9,12,15,27,33]
    rect = patches.Rectangle((wind_speed_to_reynolds(6.5),0.6), wind_speed_to_reynolds(20-6.5), 1.4-0.6, linewidth=1, edgecolor='none', facecolor='black', alpha=0.3)
    ax3.add_patch(rect)
    for i,inc in enumerate(increments):
        df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel\\'+sample+'\\'+sample+'_'+str(inc)+'_AVG.csv')
        Re = df['Re']
        Cd = df['Cd']
        ax3.plot(Re,Cd,marker='o',color=cm(inc/increments[-1]),markersize=markersize,linewidth=linewidth)
    # Define the normalization for the color mapper
    norm = Normalize(vmin=0, vmax=increments[-1]*1e-1/0.33)
    cbar2 = fig.add_cm_plot(ax3,cm,norm,ticks=np.arange(0,12,2),label='Strain, $\\epsilon$ (%)')
    ax_list.append(ax3)
    # ax3.legend(fontsize=5)

    # Generate the grid spacing figure
    # fig.figure_grid()



    # Show the main figure with all the axes
    for j,ax in enumerate(ax_list):
        # fig.set_ticks(ax,None,[0.6,1.4])

        # Create a secondary x-axis
        secax = ax.secondary_xaxis('top', functions=(reynolds_to_wind_speed, wind_speed_to_reynolds))
        fig.set_ticks(secax,[10,20,30],None)
        # Turn off x-axis ticks
        
        if j == 0:
            fig.set_ticks(ax,[5e4,1E5],[0.8,1,1.2])
            secax.set_xlabel('Wind Speed, $U$ (m/s)',fontsize=fig.fontsize_label,fontproperties=fig.font,labelpad=5)
            ax.tick_params(axis='x', which='both', labelbottom=False)
            secax.tick_params(axis='x', which='both', labeltop=True,pad=3)
            fig.set_ticks(secax,[10,20,30],None)
        elif j == 1:
            secax.tick_params(axis='x', which='both', labeltop=False,pad=3,labelsize=fig.fontsize_tick)
            # ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            # ax.xaxis.get_offset_text().set_position((1.3, 0))  # (x, y) offset in axes coordinates
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((0, 0))  # Force scientific notation
            ax.xaxis.set_major_formatter(FuncFormatter(scientific))
            # Set the font size of the scientific notation (offset text)
            ax.xaxis.get_offset_text().set_fontsize(fig.fontsize_tick)
            fig.set_ticks(ax,[5e4,1E5],[0.8,1,1.2])
            ax.set_xlabel('Reynolds Number, $Re$',fontsize=fig.fontsize_label,fontproperties=fig.font,labelpad=5)
        else:
            secax.tick_params(axis='x', which='both', labeltop=False,pad=3)
            fig.set_ticks(ax,[5e4,1E5],[0.8,1,1.2])
            ax.tick_params(axis='x', which='both', labelbottom=False)
        ax.set_xlim([1*10**4, 0.12*10**6])
        ax.set_ylim([0.6, 1.4])
        # Set the fontsize of the colorbar's ticks

        ax.spines['right'].set_visible(True)
        
        fig.set_tick_default(ax)
        fig.set_tick_default(secax,pad=1)

    cbar1.ax.tick_params(labelsize=fig.fontsize_tick)
    cbar2.ax.tick_params(labelsize=fig.fontsize_tick)
    cbar1.set_ticks([0,5,10,15])
    cbar2.set_ticks([0,5,10])
    cbar1.set_label('Strain, $\\epsilon$ (%)', labelpad=3)  # Set label with padding = 0
    cbar2.set_label('Strain, $\\epsilon$ (%)', labelpad=3)  # Set label with padding = 0
    # Process each image with a different color
    color_list = ['#333333','#EC008C']
    for idx, img_name in enumerate(['KNT','AT24']):
        size = 0.18
        label = fig.fig.add_axes([0.325, 0.34-(idx)*0.28, size, size])
        fig.circular_label_image(path+'figure_3_c\\'+img_name,label,color_list[idx],rotate=False,padding_factor=0.06,linewidth=2.25,file_type='.png')       


        

    # Add indexing to graph
    fig.add_index(padx=4,pady=-2) # Indexing in pts    # Saving of figure
    fig.savefig(path+'figure_3.pdf',dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_3()