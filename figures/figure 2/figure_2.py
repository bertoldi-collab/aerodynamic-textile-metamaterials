import sys
import os
import numpy as np
import pandas as pd
# Get the path of the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent of the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, "../"))
grandparent_dir = os.path.abspath(os.path.join(current_dir, "../../"))

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)
sys.path.append(grandparent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt
from dfplt import aspect_convert
# Import the function from the specified file
from data.scanning_cylinderical.post_process import load_heatmap

def figure_2():
    figsize = [542.03,401.24]
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    fig.fontsize_index = 12

    path = os.getcwd()

    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)

    # Figure 2 (a): 
    figsize_a = [110,200.62]
    rect_a = aspect_convert(0,figsize_a[1]/figsize[1],figsize_a,figsize)
    ax1 = fig.plot_pdf(rect_a,path+"\\figures\\figure 2\\figure_2_a\\PDF\\figure_2_a.pdf")

    # Figure 2 (b): 
    figsize_b = [162.61,200.62]
    rect_b = aspect_convert(0.15,figsize_b[1]/figsize[1],figsize_b,figsize)
    ax2 = fig.plot_pdf(rect_b,path+"\\figures\\figure 2\\figure_2_b\\PDF\\figure_2_b.pdf")

    # # Figure 2 (c)
    ssx = 0.06
    ssy = 0.09
    sample_list = ['AT10','AT9','AT24','AT11','AT15']
    fig.subplots.append([0.3,0.5,0.66,0.5])
    itr_list = [12,1]
    ax_d = []
    # Add circular images
    for idx, img_name in enumerate(sample_list):
        size = 0.1
        label = fig.fig.add_axes([(idx/5)*0.55+0.43, 0.894, size, size])
        fig.circular_label_image(path+'\\figures\\figure 2\\figure_2_c\\'+img_name,label,fig.colors[idx],rotate=False,padding_factor=0.03,linewidth=2.25)  
    # Add heatmaps     
    xlim_list =[[25,55]]
    ylim_list = [[0,30]]
    for row in [0,1]:
        for col in [0,1,2,3,4]:
            sample_name = sample_list[col]
            folder_name = 'data\\scanning_cylinderical\\'+sample_name+'\\'
            rect = [0.35+col*0.11, 0.56+0.15*(row), 0.222, 0.222]
            
            ax_i = fig.plot([rect[0]+ssx/2+0.02, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel=None, ylabel=None,pad=3)
            if sample_name != 'AT13':
                indx = itr_list[row]
                spacing = np.linspace(0,33,int(33/3+1))
            else:
                spacing = np.linspace(0,34,int(34/2+1))
                indx = 17-row*16
                spacing = np.delete(spacing,9)

            load_heatmap(ax_i,sample_name,indx,folder_name,spacing,regenerate=False)
            ax_i.set_xlim([25,55])
            ax_i.set_ylim([0,30])
    ax_c = fig.plot([0+ssx/2, 0+0.05, 0.6-ssx, 0.5-0.15])
    data = np.full((10, 10), np.nan)
    fig.plot_heatmap(data,ax_c,[0,0.5,1,1.5],'Surface Height, $z$ (mm)',pos=[0.15+0.43, 0.57, 0.3, 0.02],NaN_colorbar=True,pad_cb=0.025,vmax =1.75,vmin=0)
    # ax.set_xlim((30,60))
    # ax.set_ylim((0,30))
    # Add in Scale bar
    ax3_2 = fig.plot([0.427, 0.58, 0.11, 0.02])
    ax3_2.plot([0,24],[0,0],color='k',linewidth=1)
    ax3_2.plot([0,0],[-2,2],color='k',linewidth=0.9)
    ax3_2.plot([24,24],[-2,2],color='k',linewidth=1)
    ax3_2.set_xlim([-3.75,78.75])
    ax3_2.axis('off')
    # ax3_2.set_xlabel('20mm', fontsize=fig.fontsize_label,fontproperties=fig.font)
    fig.fig.text(0.421, 0.57, '10mm', fontsize=fig.fontsize_label, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
    


    # Figure 2 (e)
    ssx = 0.05
    ssy = 0.09
    rect = [0.5, 0.025, 0.4, 0.5]
    fig.subplots.append(rect)
    ax3 = fig.plot([rect[0]+ssx/2+0.04, rect[1]+ssy/2, rect[2]-ssx, rect[3]-ssy], xlabel="Strain, $\epsilon$ (%)", ylabel="Surface Roughness, $\overline{k}$ (mm)",pad=3)
    sample_list = ['AT10','AT9','AT24','AT11','AT15']
    # fig = ['#ec008c','#ec008c','#264653','#264653','#2A9D8F','#2A9D8F','#E9C46A','#E9C46A','#F4A261','#F4A261']
    for i,sample in enumerate(sample_list):
        df = fig.load_csv(path+'\\data\\scanning_cylinderical\\'+sample+'\\roughness_displacement.csv')
        inc = df['Displacement']
        depths = df['Roughness']
        ax3.plot(inc/3.3,depths,marker='o',color=fig.colors[i])
        D = 57.15

        df = fig.load_csv(path+'\\data\\scanning_cylinderical\\'+sample+'\\{}_SD.csv'.format(sample))
        inc = df['Displacement']
        roughness_SD = df['Roughness_SD']
        ax3.errorbar(inc/3.3, depths, yerr=roughness_SD, fmt='none', ecolor=fig.colors[i], capsize=0,alpha=0.5)
    fig.set_ticks(ax3,[0,5,10],[0.2,0.4,0.6,0.8,1,1.2])
    # Create the secondary y-axis and scale the y-values by dividing by D
    ax3_twin = ax3.twinx()  # instantiate a second axes that shares the same x-axis
    ax3_twin.set_ylim(ax3.get_ylim()[0] / D, ax3.get_ylim()[1] / D)  # Scale y-limits by D
    ax3_twin.set_ylabel(r'Roughness Coefficient, $\overline{k}/D$',fontproperties=fig.font,fontsize=fig.fontsize_label,labelpad=5)  # Set the secondary y-axis label
    for label in ax3_twin.get_yticklabels():
        label.set_fontproperties(fig.font)
    ax3_twin.tick_params(axis='y', labelsize=fig.fontsize_tick,direction='in')
    ax3_twin.set_yticks([0.01,0.02])
    ax3_twin.spines['top'].set_visible(False)

    # Add indexing to graph
    fig.add_index(padx=1,pady=-1) # Indexing in pts

    # Saving of figure
    fig.savefig(path+'\\figures\\figure 2\\figure_2_temp.pdf',dpi=500)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_2()