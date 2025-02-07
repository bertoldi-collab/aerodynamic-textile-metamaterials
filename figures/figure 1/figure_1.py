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

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt
from dfplt import aspect_convert

def figure_1():

    figsize = [510.236,377.704]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    # OS Cwd
    path = os.getcwd()+"\\figures\\figure 1\\"
    parentfolder = os.getcwd()

    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)
    fig.fontsize_index = 12
    linewidth_graphics = 1

    # Figure 1 (a): 
    figsize_a = [311.56,188.852]
    rect_a = aspect_convert(0,1-figsize_a[1]/figsize[1],figsize_a,figsize)
    ax1 = fig.plot_pdf(rect_a,path+"figure_1_a\\PDF\\figure_1_a.pdf") 

    # Figure 1 (b)
    figsize_b = [178.6762,188.852]
    rect_b = aspect_convert(1-figsize_b[0]/figsize[0],1-figsize_b[1]/figsize[1],figsize_b,figsize)
    ax2 = fig.plot_pdf(rect_b,path+"figure_1_b\\PDF\\figure_1_b.pdf") 
    fig.fig.text(0.69, 0.53, '$\\epsilon = 0\\%$', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
    fig.fig.text(0.8, 0.53, '$\\epsilon = 10\\%$', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)

    # Add in Scale bar
    ax2_2 = fig.plot([0.675, 0.89, 0.1, 0.02])
    ax2_2.plot([0,0.3],[0,0],color='k',linewidth=linewidth_graphics)
    ax2_2.plot([0,0],[-2,2],color='k',linewidth=linewidth_graphics)
    ax2_2.plot([0.3,0.3],[-2,2],color='k',linewidth=linewidth_graphics)
    ax2_2.set_xlim([-0.01,1])
    ax2_2.axis('off')
    # ax3_2.set_xlabel('20mm', fontsize=fig.fontsize_label,fontproperties=fig.font)
    fig.fig.text(0.665, 0.94, '20mm', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)

    # Plot your heat-map
    rect = [0.91, 0.53, 0.15, 0.24]
    ax_c = fig.plot(rect)
    data = np.full((10, 10), np.nan)
    cbar = fig.plot_heatmap(data,ax_c,[-4,-2,0,2,4],'Surface Height, $z$ (mm)',pos=rect,NaN_colorbar=True,pad_cb=3,vmax =4.36556,vmin=-4.36556,vertical=True)
    cbar.ax.tick_params(width=0.7)

    # Format the colorbar ticks
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Figure 1 (c)
    ssx = 0.07
    ssy = 0.00
    fig.subplots.append([0, 0, 0.6, 0.5])
    ax3 = fig.plot([0+ssx/2, 0+0.05, 0.6-ssx, 0.5-0.15])
    sample_list = ['AT10','AT9','AT24','AT11','AT15']
    index_offset = [1,3,1,2,1]
    data = None
    df_min = 1E20
    # File name
    filename = 'data.pickle'
    # Check if file exists
    if os.path.exists(path+'\\figure_1_c\\'+filename):
        # Load the data from the pickle file
        with open(path+'\\figure_1_c\\'+filename, 'rb') as file:
            data,df_min = pickle.load(file)
        print("Loaded data from data.pickle")
    else:
        for i,sample in enumerate(sample_list):
            df = fig.load_csv(path+'\\figure_1_c\\'+sample+'_depth_{}.csv'.format(15+index_offset[i]))
            data = pd.concat([data, df], axis=1)

            if df_min > np.min(np.shape(df)[0]):
                df_min = np.min(np.shape(df)[0])
        # Save the data to the pickle file
        with open(path+'\\figure_1_c\\'+filename, 'wb') as file:
            pickle.dump((data,df_min), file)
        print("Data pickled and saved to data.pickle")
    # Slice both DataFrames to the minimum number of rows
    data = data.iloc[:df_min, :]
    # Plot your heat-map
    fig.plot_heatmap(data,ax3,[-4,-2,0,2,4],'Surface Height, $z$ (mm)',pos=[0.15, 0.07, 0.3, 0.02],NaN_colorbar=True,pad_cb=0.025,vmax =4.36556,vmin=-4.36556)
    # Process each image with a different color
    for idx, img_name in enumerate(sample_list):
        size = 0.1
        label = fig.fig.add_axes([(idx/5)*0.528+0.041, 0.41, size, size])
        fig.circular_label_image(path+'figure_1_c\\'+img_name,label,fig.colors[idx],rotate=False,padding_factor=0.03,linewidth=2.25)       

    # Add in Scale bar
    ax3_2 = fig.plot([0.035, 0.08, 0.11, 0.02])
    ax3_2.plot([0,20],[0,0],color='k',linewidth=linewidth_graphics)
    ax3_2.plot([0,0],[-2,2],color='k',linewidth=linewidth_graphics)
    ax3_2.plot([20,20],[-2,2],color='k',linewidth=linewidth_graphics)
    ax3_2.set_xlim([-3.75,78.75])
    ax3_2.axis('off')
    # ax3_2.set_xlabel('20mm', fontsize=fig.fontsize_label,fontproperties=fig.font)
    fig.fig.text(0.03, 0.06, '20mm', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
        
    # Add in 10% Strain
    yarr = 90.909091
    ax3_3 = fig.plot([0.035, 0.052, 0.14, 0.404])
    ax3_3.plot([0,20],[yarr,yarr],':',color='k',linewidth=linewidth_graphics)
    ax3_3.plot([0,20],[100,100],':',color='k',linewidth=linewidth_graphics)
    ax3_3.plot([3,3],[100,110],'-',color='k',linewidth=linewidth_graphics)
    ax3_3.plot([3,3],[yarr-10,yarr],'-',color='k',linewidth=linewidth_graphics)
    # ax3_3.plot([3,3],[5.8+3,5.8],'-',color='k',linewidth=linewidth_graphics)
    plt.annotate('', xy=(3, 100), xytext=(3, 101),
             arrowprops=dict(facecolor='black',headwidth=4,headlength=4,lw=0.5))
    # ax3_3.plot([3,3],[-3,0],'-',color='k',linewidth=linewidth_graphics)
    plt.annotate('', xy=(3, yarr), xytext=(3, yarr-1),
            arrowprops=dict(facecolor='black',headwidth=4,headlength=4,lw=0.5))

    ax3_3.set_xlim([-3.75,78.75])
    ax3_3.set_ylim([-20,120])
    ax3_3.axis('off')
    # ax3_2.set_xlabel('20mm', fontsize=fig.fontsize_label,fontproperties=fig.font)
    fig.fig.text(0.06, 0.35, '$\\epsilon = 10\\%$', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
    

    
    # Figure 1 (d)
    ssx = 0.08
    ssy = 0.10
    rect = [0.59, 0.01, 0.4, 0.5]
    fig.subplots.append([0.58, 0.0, 0.4, 0.5])
    ax4 = fig.plot([rect[0]+ssx/2+0.03, rect[1]+ssy/2+0.01, rect[2]-ssx, rect[3]-ssy], xlabel="Strain, $\epsilon$ (%)", ylabel="Surface Roughness, $\overline{k}$ (mm)",pad=3)
    sample_list = ['AT10','AT9','AT24','AT11','AT15']
    # Index Offset: Due to planar nature of membrane samples and boundary condition affects, it is possible for the sample to 'begin' slack 
    # at initial strain values. Index offset will re-zero the strain values before the initial rise of results
    index_offset = [1,3,1,2,1]
    for i,sample in enumerate(sample_list):
        ind = index_offset[i]
        df = fig.load_csv(os.getcwd+'\\data\\scanning_planar\\'+sample+'_PL\\'+sample+'.csv')
        inc = df['increments']
        depths = df['k_average']
        ax4.plot((inc[ind:ind+16]-inc[ind])/1.5,depths[ind:ind+16],marker='o',color=fig.colors[i])

        # Standard Deviation Values
        df = fig.load_csv(path+'figure_1_d\\'+'{}_SD.csv'.format(sample))
        inc = df['increments']
        roughness_SD = df['k_average']
        ax4.errorbar((inc[ind:ind+16]-inc[ind])/1.5,depths[ind:ind+16], yerr=roughness_SD[ind:ind+16], fmt='none', ecolor=fig.colors[i], capsize=0,alpha=0.5)
    fig.set_ticks(ax4,[0,5,10],[1,3,5,7])

    coordinates =[[0.85,0.102],[0.74,0.17],[0.92,0.1],[0.87,0.26],[0.87,0.35]]
    # Process each image with a different color
    for idx, (coord,img_name) in enumerate(zip(coordinates,sample_list)):
        size = 0.06
        label = fig.fig.add_axes([coord[0], coord[1], size, size])
        fig.circular_label_image(path+'figure_1_c\\'+img_name,label,fig.colors[idx],rotate=False,padding_factor=0.05,linewidth=2)       


    # Figure 1 (d_ii): 
    figsize_dii = [39.5,54.598]
    rect_d_ii = aspect_convert(0.675,0.31,figsize_dii,figsize)
    axd_ii = fig.plot_pdf(rect_d_ii,path+"figure_1_d\\PDF\\figure_d_ii.pdf",add_index=False) 

    # Generate the grid spacing figure
    fig.figure_grid()

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # # Show the main figure with all the axes
    # fig.show()

    # Saving of figure
    fig.savefig(path+'figure_1.pdf',dpi=800)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_1()