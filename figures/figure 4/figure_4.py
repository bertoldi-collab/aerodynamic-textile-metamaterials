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
import matplotlib.patches as patches
from shapely.geometry import MultiPoint, Polygon
from shapely.ops import cascaded_union
import alphashape
from matplotlib.collections import PolyCollection
from scipy.spatial import Delaunay
from scipy.stats import zscore

# Get the path of the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent of the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, "../"))

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt
from dfplt import aspect_convert

def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n, 2) points
    :param alpha: alpha parameter to control the shape
    :return: set of simplices defining the alpha shape
    """
    if len(points) < 4:
        # If the number of points is less than 4, can't construct an alpha shape
        return []

    tri = Delaunay(points)
    valid_triangles = []
    for ia, ib, ic in tri.simplices:
        pa, pb, pc = points[ia], points[ib], points[ic]
        # Compute circumradius
        a = np.linalg.norm(pa - pb)
        b = np.linalg.norm(pb - pc)
        c = np.linalg.norm(pc - pa)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = (a * b * c) / (4.0 * area)
        if circum_r < 1.0 / alpha:
            valid_triangles.append([ia, ib, ic])

    return valid_triangles

def detect_outliers(points, threshold=3):
    """
    Detect and remove outliers from a set of 2D points using the z-score method.
    :param points: np.array of shape (n, 2) points
    :param threshold: z-score threshold to classify points as outliers
    :return: points without outliers
    """
    z_scores = np.abs(zscore(points, axis=0))
    mask = (z_scores < threshold).all(axis=1)  # Retain points with all z-scores below threshold
    return points[mask]

def figure_4(dpi=500):
    def scientific(x, pos):
        return f'{x:.0e}'

    def reynolds_to_wind_speed(x):
        return x*dynamic_viscosity/diameter_cylinder

    def wind_speed_to_reynolds(x):
        return x*diameter_cylinder/dynamic_viscosity
    figsize = [510.236,360]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    # colors =['#CDCDCD','#333333']
    # cm = fig.custom_cm(colors)
    cm = plt.get_cmap('jet')

    # Constant:
    diameter_cylinder = 0.05715     # 2.25inch cylinder (mm)
    dynamic_viscosity = 1.516*10**-5 # m^2/s  (20 degree C)

    # OS Cwd
    path = os.getcwd()+"\\figures\\figure 4\\"
    ax_list = []
    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    
    # figure 4 (a)
    figsize_a = [100,216.81]
    rect_a = aspect_convert(0,1-figsize_a[1]/figsize[1],figsize_a,figsize)
    ax1 = fig.plot_pdf(rect_a,path+"figure_4_a\\PDF\\figure_4_a.pdf") 

    # figure 4 (b)
    ssx = 0.13
    ssy = 0.13
    rect = [0.18, 0.5, 0.45, 0.5]
    
    # for width in  np.linspace(10,21,12):
	# 	for height in np.linspace(12,43,32):
	# 		for per_alpha in np.linspace(0.0,0.35,15):
	# 			alpha = np.arctan((height/2)/(width*per_alpha))
	# 			tasks.append([width,height,alpha,cyl_offset,disp_load])
    fig.subplots.append([0.16,0.6,0.31,0.4])
    ax6 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel="Area, $w \\times h$ ($\mathrm{mm}^2$)", ylabel="Max Roughness, $k_{max}$ (mm)",pad=3)
    df = fig.load_csv(path+'figure_4_a\\sampling_data.csv')
    failed_list = []
    sample_list = ['AT12','AT13','AT14','AT16','AT17','AT18','AT19','AT20','AT21','AT22','AT23','AT24','AT25']
    cc_excluded = 0
    colors = ['#ec008c','#F4A261','#2A9D8F','#264653','#E9C46A']
    n1_list = []
    n2_list = []
    alpha_list = []
    collapsed = 0
    for i,sample in enumerate(np.array(df)):
        # print(sample)
        folder_name = sample[0]
        area = sample[1]
        width = sample[2]
        height = sample[3]
        alpha = sample[4]
        x_offset = np.round((height / 2) / (width * np.tan(alpha*np.pi/180)),4)
        ep = sample[5]
        k = sample [6]
        N = sample[7]
        if  x_offset >= 0.10:
            alpha_list.append(alpha)
            # if alpha < 1.45 and height > 16 and width < 16:
            if folder_name in ['AT24','AT23','AT19','AT22','AT18']:
                index = np.where(folder_name==np.array(['AT24','AT23','AT19','AT22','AT18']))[0][0]
                if N ==1:
                    ax6.scatter(area, k, alpha=1,color =colors[index],marker='D',zorder=3,s=60)
                else:
                    ax6.scatter(area, k, alpha=1,color =colors[index],marker='s',zorder=3,s=60)
                
                # Label the point with the number '1'
                ax6.text(
                    area - 70, k+0.08, 'G'+str(index+1),  # Adjust the position with offsets
                    fontsize=fig.fontsize_index,            # Set font size
                    color=colors[index],            # Set font color
                    fontweight='bold'       # Optional: Make it bold
                )
            elif ep != 0.1:
                if ep < 0.1:
                    failed_list.append([folder_name])
                else:
                    cc_excluded = cc_excluded + 1
            else:
                if N == 1:
                    ax6.scatter(area, k, alpha=0.12,color ='#77CDFF',marker='D',zorder=2)
                    n1_list.append([area,k])
                else:
                    ax6.scatter(area, k, alpha=0.05,color ='#606060',marker='s',zorder=1)
                    n2_list.append([area,k])
                    collapsed = collapsed + 1
        else:
            cc_excluded = cc_excluded + 1
    
    ax6.set_xlim([30,960])
    ax6.set_ylim([0.1,2.2])
    fig.set_ticks(ax6,[100,300,500,700,900],[0.5,1.0,1.5,2.0],length=3)

    # axb_cm = fig.plot([0.3, 0.8, 0.15, 0.02], xlabel=None, ylabel=None,pad=0)
    # fig.plot_heatmap(data,axb_cm,[-4,-2,0,2,4],' (mm)',pos=rect,NaN_colorbar=True,pad_cb=3,vmax =,vmin=0,vertical=True)

    print('Number failed: {}'.format(len(failed_list)))
    print('Percentage Failed: {}'.format(len(failed_list)/(i-len(sample_list)-cc_excluded)))
    print('Percentage Collapsed: {}'.format(collapsed/(i-len(sample_list)-cc_excluded)))
    print('Total Simulations: {}'.format(i-len(sample_list)-cc_excluded))
    print('Alpha Min/Max: {},{}'.format(np.min(np.array(df)[:,4]),np.max(np.array(df)[:,4])))
    print('Alpha Min: {}'.format(i-len(sample_list)-cc_excluded))
    

    # figure 4 (c)
    rect = [0.58, 0.5, 0.4, 0.5]
    fig.subplots.append([0.56, 0.6, 0.4, 0.4])
    ax4 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel="Strain, $\\epsilon$ (%)", ylabel="Tared Roughness, $\overline{k}-\overline{k}_0$ (mm)",pad=3)
    sample_list = ['AT24','AT23','AT19','AT22']
    colors = ['#ec008c','#2A9D8F','#F4A261','#264653','#E9C46A']
    baseline=0
    # Create custom lines for the legend
    experiment_line = plt.Line2D([0], [0], color='black', marker='o', label='Experiment')
    simulation_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Simulation')
    markersize = 2.5
    linewidth = 1.5
    baseline_list = []
    for i,sample in enumerate(sample_list):
        df_exp = fig.load_csv(os.getcwd()+'\\data\\scanning_cylinderical\\'+sample+'\\'+'roughness_displacement.csv')
        strain_exp = df_exp.iloc[:, 0]
        roughness_exp = df_exp.iloc[:, 1]

        df_sim = fig.load_csv(os.getcwd()+'\\data\\FEM_cylinderical_PBC\\'+sample+'\\'+'depth.csv')
        strain_sim = df_sim.iloc[:, 0]
        roughness_sim = df_sim.iloc[:, 1]

        # Plot Experimental data
        baseline= roughness_exp[0]
        baseline_list.append(roughness_exp[0])
        ax4.plot(strain_exp/3.3,roughness_exp-baseline,marker='o',color=colors[i],markersize=markersize)
        # Add error bars separately using ax.errorbar
        df_sd = fig.load_csv(os.getcwd()+'\\data\\scanning_cylinderical\\'+sample+'\\'+sample+'_SD.csv')
        roughness_SD = df_sd.iloc[:, 1]
        ax4.errorbar(strain_exp/3.3, roughness_exp-baseline, yerr=roughness_SD[0:len(roughness_exp)], fmt='none', ecolor=colors[i], capsize=0,alpha=0.5)
        
        # Plot Simulation data
        ax4.plot(strain_sim*100,roughness_sim,'--',color=colors[i])


    print(baseline_list)
    print(np.mean(baseline))


    # Add the custom legend to the plot
    markersize = 2.5
    linewidth = 1.5
    D = 57.15
    ax4.legend(handles=[experiment_line, simulation_line],fontsize=fig.fontsize_tick-1,loc='lower right')

    ax4.set_ylim([-0.1,2.4])
    fig.set_ticks(ax4,[0,5,10],[0,1,2],length = 3)
    ax3_twin = ax4.twinx()  # instantiate a second axes that shares the same x-axis
    ax3_twin.set_ylim(ax4.get_ylim()[0] / D, ax4.get_ylim()[1] / D)  # Scale y-limits by D
    ax3_twin.set_ylabel(r'Roughness Coefficient, $k/D$',fontproperties=fig.font,fontsize=fig.fontsize_label,labelpad=5)  # Set the secondary y-axis label
    for label in ax3_twin.get_yticklabels():
        label.set_fontproperties(fig.font)
    ax3_twin.tick_params(axis='y', labelsize=fig.fontsize_tick,direction='in',length =3)
    ax3_twin.set_yticks([0.01,0.02,0.03,0.04])
    ax3_twin.spines['top'].set_visible(False)
    data = np.full((10, 10), np.nan)


    # figure 4 (c)
    ssx = 0.13
    ssy = 0.13
    baseline=0
    # Create custom lines for the legend
    experiment_line = plt.Line2D([0], [0], color='black', marker='o', label='Experiment')
    simulation_line = plt.Line2D([0], [0], color='black', linestyle='--', label='Simulation')
    sample_list = ['AT24','AT23','AT19','AT22']
    increments = [0,3,9,12,15,27,33]
    fig.subplots.append([0,0.0,1,0.5])
    color_list = [['#fbcce8','#EC008C'],['#fdecdf','#F4A261'],['#d4ebe9','#2A9D8F'],['#d4dadd','#264653']]
    df_smooth_cylinder = fig.load_csv(path+'smooth_cylinder.csv')
    Re_sc = df_smooth_cylinder['Re']
    Cd_sc = df_smooth_cylinder['Cd']
    for i, sample in enumerate(sample_list):
        rect = [0.02+i*0.21, 0.02, 0.32, 0.44]
        if i == len(sample_list) - 1:
            rect[2] = 0.36
        ax_i = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel=None, ylabel=None, pad=3)
        if i == 0:
            ax_i.set_ylabel("Drag Coefficient $C_d$", fontsize=fig.fontsize_label, fontproperties=fig.font, labelpad=3)
        fig.set_ticks(ax_i, [5E4, 1E5], [0.6, 0.8, 1, 1.2], length=3)
        ax_i.tick_params(axis='y', which='both', labelleft=(i == 0), length=3)
        cm = plt.get_cmap('jet')
        for inc in increments:
            df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel\\'+sample+'\\{}_{}_AVG.csv'.format(sample, int(inc)))
            Re = df['Re']
            Cd = df['Cd']
            ax_i.plot(Re, Cd, marker='o', color=cm(inc/increments[-1]), markersize=markersize)
        ax_i.set_xlim([1*10**4, 0.12*10**6])
        ax_i.set_ylim([0.6, 1.4])
        ax_i.set_xticks([5e4, 1e5])
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 0))  # Force scientific notation
        ax_i.xaxis.set_major_formatter(FuncFormatter(scientific))
        rect = patches.Rectangle((wind_speed_to_reynolds(6.5), 0.6), wind_speed_to_reynolds(20-6.5), 1.4-0.6, linewidth=1, edgecolor='none', facecolor='black', alpha=0.3)
        ax_i.add_patch(rect)
        ax_i.plot(Re_sc, Cd_sc, '--k', linewidth=linewidth, label='Smooth Cylinder [7]',zorder=0)
        secax = ax_i.secondary_xaxis('top', functions=(reynolds_to_wind_speed, wind_speed_to_reynolds))
        fig.set_ticks(secax,[10,20,30],None,length = 3)
        secax.set_xlabel(None,fontsize=fig.fontsize_label,fontproperties=fig.font,labelpad=3)
        secax.tick_params(axis='x', which='both', labeltop=True,pad=3,length =3)

        
    # Add common x-axis label
    fig.fig.text(0.5, 0.02, 'Reynolds Number, ${Re}$', ha='center', fontsize=fig.fontsize_label, fontproperties=fig.font)
    fig.fig.text(0.5, 0.455, 'Wind Speed, $U$ (m/s)', ha='center', fontsize=fig.fontsize_label, fontproperties=fig.font)

    norm = Normalize(vmin=0, vmax=increments[-1]*1e-1/0.33)
    fig.add_cm_plot(ax_i, cm, norm, ticks=np.arange(0, 12, 2), label='Strain, $\\epsilon$ (%)')

    for idx, img_name in enumerate(sample_list):
        size = 0.06
        label = fig.fig.add_axes([(idx/4)*0.835+0.13, 0.13, size, size])
        fig.circular_label_image(path+'figure_4_c\\'+img_name, label, color_list[idx][1], rotate=False, padding_factor=0.03, linewidth=1.85)


    # # Generate the grid spacing figure
    # fig.figure_grid()

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1.2) # Indexing in pts    # Saving of figure
    fig.savefig(path+'figure_4_temp.pdf',dpi=1000)

    

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_4()