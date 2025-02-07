import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

# Get the path of the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent of the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, "../"))

# Add the grandparent directory to sys.path
sys.path.append(parent_dir)

# Now you can import dfplt from dfplt.py located in the grandparent directory
from dfplt import dfplt


def figure_5():
    def reynolds_to_wind_speed(x):
        return x*dynamic_viscosity/diameter_cylinder

    figsize = [249.5,460]

    # Create a large main figure
    fig = dfplt(fignum=2,figsize=figsize) # Size in pt's
    x_spacing =0.15
    ssx = 0.0
    ssy = 0.0
    spacing = 0.025
    path = os.getcwd()+"\\figures\\figure 5\\"
 
    # Constant:
    diameter_cylinder = 0.05715     # 2.25inch cylinder (mm)
    dynamic_viscosity = 1.516*10**-5 # m^2/s  (20 degree C)

    markersize = 4
    linewidth = 2

    # Figure 5 (a)
    ssx = 0.09
    ssy = 0.09
    rect =[0.15, 0.66, 0.85, 0.33]
    

    ax1 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2, rect[2]-ssx-0.1, rect[3]-ssy], xlabel='Wind-speed, $U$ (m/s)', ylabel='Min. Drag Coefficient, $C_{d,min}$',pad=3)
    colors = ['#7E7E7E','#01B1CE','#F4A261','#EC008C']
    line5 = plt.Line2D([0], [0], color=colors[1], linestyle='-',linewidth=4, label='$\\epsilon = 5\\%$')
    line10 = plt.Line2D([0], [0], color=colors[2], linestyle='-',linewidth=4, label='$\\epsilon = 10\\%$')
    dynamic_line = plt.Line2D([0], [0], color=colors[0], linestyle='-',linewidth=4, label='$\\epsilon = 0\\%$')
 

    rect = [0.15, 0.39, 0.85, 0.33]
    fig.subplots.append([0.03, 0.66, 0.85, 0.33])
    sample_list = ['AT23']
    colors = ['k','#01B1CE','#F4A261','#EC008C']
    increments =[0,15,33]
    markersize = 5
    linewidth = 2.5
    n = 14
    speeds = np.linspace(0,n-1,n,dtype=int)
    cm = plt.get_cmap('jet')
    for j,sample in enumerate(sample_list):
        e_optimal = 1E9*np.ones((n,3))
        for speed in speeds:
            for i,inc in enumerate(increments):
                df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel\\'+sample+'\\'+sample+'_'+str(inc)+'_AVG.csv')
                Re = df['Re']
                Cd = df['Cd']
                ax1.plot(reynolds_to_wind_speed(Re),Cd,marker='o',color=colors[i],markersize=markersize, alpha = 0.02)
                if inc == 15:
                    ax1.plot(reynolds_to_wind_speed(Re[7:14]),Cd[7:14],color=colors[1],marker='o',markersize=markersize,linewidth=linewidth,alpha=1,zorder =0)
                elif inc == 33:
                    ax1.plot(reynolds_to_wind_speed(Re[0:8]),Cd[0:8],color=colors[2],marker='o',markersize=markersize,linewidth=linewidth,alpha=1)
                if Cd[speed] < e_optimal[speed,0]:
                    e_optimal[speed,0] = Cd[speed]
                    e_optimal[speed,1] = inc
                    e_optimal[speed,2] = Re[speed]

    ax1.set_ylim([0.6,1.25])
    ax1.set_xlim([6,20.5])
    ax1.legend(loc='lower right',handles=[line5, line10,dynamic_line],fontsize=fig.fontsize_label,frameon=False)

    fig.set_ticks(ax1,[10,15,20],[0.8,1.0,1.2],length = 3)
    rect = patches.Rectangle((6.5,0.6), 20-6.5, 1.4-0.6, linewidth=1, edgecolor='none', facecolor='black', alpha=0.12)
    ax1.add_patch(rect)
    for idx, img_name in enumerate(sample_list):
        size = 0.15
        label = fig.fig.add_axes([(idx/4)*0.835+0.3, 0.69, size, size])
        fig.circular_label_image(path+'figure_5_a\\'+img_name,label,color='#2A9D8F',rotate=False,padding_factor=0.03,linewidth=2)       
    
    # Figure 5 b  - Dynamic
    ssx = 0.09
    ssy = 0.09
    rect = [0.15, 0.34, 0.85, 0.33]
    fig.subplots.append([0.03, 0.345, 0.85, 0.33])
    ax7 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2+0.01, rect[2]-ssx-0.1, rect[3]-ssy], xlabel="Time, $t$ (min)", ylabel="Wind-Speed, $V$ (m/s)",pad=3)
    # Dummy plots for legend (with 'None' values)
    ax7.plot([], [], label="5% Strain", color='#01B1CE',markersize=1.5)
    ax7.plot([], [], label="10% Strain", color='#F4A261',markersize =1.5)
    sample = 'profile_velocities'
    df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel_active\\'+sample+'.csv')
    time = df['time']
    velocity = df['velocity']
    ax7.plot(time/60,velocity,color='k',linewidth=linewidth)
    fig.set_ticks(ax7,[0,5,10,15],[10,15])
    ax7.set_xlim([0,10])
    ax7.set_ylim([4,18])
    rect = patches.Rectangle((0,6.5), 10, 13.5, linewidth=1, edgecolor='none', facecolor='black', alpha=0.2)
    ax7.tick_params(axis='both',direction='in', length=4)  # Set tick length to 10

    fig.fontsize_tick = 9
    fig.fontsize_label = 10
    ssx = 0
    ssy = 0
    rect = [0.35, 0.83,0.15, 0.12]
    line, = ax1.plot([14,15],[10,4.545],color='#565656',linestyle=':',linewidth=2)
    line.set_dashes([1, 1])
    sample_list = ['AT23']
    colors = ['#01B1CE','#F4A261','#EC008C']
    increments =[15,33]
    markersize = 2.5
    linewidth = 1.5
    n = 13
    speeds1 = np.linspace(0,7,8,dtype=int)
    speeds2 = np.linspace(8,n-1,n-8,dtype=int)

    # Arrows and Phrases
    fig.fig.text(0.30, 0.33+0.123, '$\\epsilon = 10\\%$\noptimum', fontsize=fig.fontsize_label-1, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
    fig.fig.text(0.64, 0.33+0.123, '$\\epsilon = 5\\%$\noptimum', fontsize=fig.fontsize_label-1, ha='left', va='top', fontproperties=fig.font)#,fontfamily=self.font_bold)
    xb = 4.7
    arrow_size = 7
    ax7.annotate('', xy=(0, xb), xytext=(0.1, xb),
            arrowprops=dict(facecolor='#F4A261',headwidth=arrow_size,headlength=arrow_size,lw=0.1,edgecolor='none'))
    ax7.annotate('', xy=(5, xb), xytext=(4.9, xb),
            arrowprops=dict(facecolor='#F4A261',headwidth=arrow_size,headlength=arrow_size,lw=0.1,edgecolor='none'))
    # ax3_3.plot([3,3],[-3,0],'-',color='k',linewidth=linewidth_graphics)
    ax7.annotate('', xy=(5, xb), xytext=(5.1, xb),
            arrowprops=dict(facecolor='#01B1CE',headwidth=arrow_size,headlength=arrow_size,lw=0.1,edgecolor='none'))
    ax7.annotate('', xy=(10, xb), xytext=(9.9, xb),
            arrowprops=dict(facecolor='#01B1CE',headwidth=arrow_size,headlength=arrow_size,lw=0.1,edgecolor='none'))
    ax7.plot([0,5],[xb,xb],'-',color='#F4A261',linewidth=1.5)
    ax7.plot([5,10],[xb,xb],'-',color='#01B1CE',linewidth=1.5)



    fig.fontsize_tick = 9
    fig.fontsize_label = 10
    # Figure 5 c  - Active Aerodynamic Test
    ssx = 0.09
    ssy = 0.09

    rect = [0.15, 0.02, 0.85, 0.33]
    fig.subplots.append([0.03, 0.02, 0.85, 0.33])
    window_size_list = [2000,2000,2000,2000]

    colors = ['k','#01B1CE','#F4A261','#EC008C']
    ax9 = fig.plot([rect[0]+ssx/2, rect[1]+ssy/2+0.01, rect[2]-ssx-0.1, rect[3]-ssy], xlabel="Time, $t$ (min)", ylabel="Force Reduction, $\Delta F_r$ (%)",pad=3)
    line5 = plt.Line2D([0], [0], color=colors[1], linestyle='-',linewidth=4, label='$\\epsilon = 5\\%$')
    line10 = plt.Line2D([0], [0], color=colors[2], linestyle='-',linewidth=4, label='$\\epsilon = 10\\%$')
    dynamic_line = plt.Line2D([0], [0], color=colors[3], linestyle='-',linewidth=4, label='$\\epsilon = 5\\leftrightarrow 10\%$')
    std_velocity_reduction_residual_list = []
    velocity_reduction_smooth_list = []
    for k,increment in enumerate([0,15,33,40]):
        velocity_list = []
        force_list = []
        force_list_std = []
    
        df = fig.load_csv(os.getcwd()+'\\data\\wind_tunnel_active\\AT23_{}_1.csv'.format(increment,i))
        window_size = window_size_list[k]

        velocity_list.append(gaussian_filter1d(df['V'], sigma=3))
        force_list.append(gaussian_filter1d(df['F'], sigma=3))
        force_list_std.append(df['F'])

        average_velocity = np.mean(velocity_list,axis=0)
        average_force = np.mean(force_list,axis=0)

        if increment == 0:
            force_0 = average_force
            velocity_0 = average_velocity
        window_size = 10000
        force_reduction = 100*(force_0-average_force)/force_0
        velocity_reduction = 100*(velocity_0-average_velocity)/velocity_0

        force_reduction_smooth = savgol_filter(np.array(force_reduction),window_size,3)
        velocity_reduction_smooth = savgol_filter(np.array(velocity_reduction),window_size,3)
        force_reduction_residual = force_reduction-force_reduction_smooth
        velocity_reduction_residual = velocity_reduction-velocity_reduction_smooth

        std_force_reduction_residual = np.std(force_reduction_residual)
        std_velocity_reduction_residual = np.std(velocity_reduction_residual)

        if increment != 0:
            plt.fill_between(time/60, force_reduction_smooth - std_force_reduction_residual, force_reduction_smooth + std_force_reduction_residual, color=colors[k], alpha=0.3)
            ax9.plot(time/60, force_reduction_smooth, color=colors[k], linewidth=linewidth, alpha=1)
            
            # Save the velocity data in lists
            velocity_reduction_smooth_list.append(velocity_reduction_smooth)
            std_velocity_reduction_residual_list.append(std_velocity_reduction_residual)

    fig.set_ticks(ax9,[0,5,10,15],[0,5,10,15,20])
    ax9.legend(loc='lower right',handles=[line5, line10,dynamic_line],fontsize=fig.fontsize_label,frameon=False)

    ax9.set_xlim(0,10)

    ax7.set_facecolor('none')
    ax9.set_facecolor('none')
    ax9.tick_params(axis='both',direction='in', length=4)  # Set tick length to 10

    # Add indexing to graph
    fig.add_index(padx=0,pady=0) # Indexing in pts    # Saving of figure
    fig.savefig(path+'figure_5.pdf',dpi=1000)

if __name__ == "__main__":
    # Basic Testing of figure generator
    figure_5()