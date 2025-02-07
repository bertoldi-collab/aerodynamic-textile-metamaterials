import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pdf2image import convert_from_path
import pandas as pd
import numpy as np
from matplotlib.cm import ScalarMappable
from PIL import Image, ImageDraw
import os

# Load the Helvetica TTF file from your local system
helvetica_path = "C:/Users/david/AppData/Local/Microsoft/Windows/Fonts/Helvetica.ttf"
helvetica_bold_path = "C:/Users/david/AppData/Local/Microsoft/Windows/Fonts/Helvetica-Bold.ttf"

# Create a font object using the TTF file
helvetica_font = fm.FontProperties(fname=helvetica_path)
helvetica_bold_font = fm.FontProperties(fname=helvetica_bold_path)

# Color Packages
colors =['#E9C46A','#2A9D8F','#ec008c','#F4A261','#264653','#01B1CE','#3FB617']
# colors = ['#01B1CE','#3FB617','#F1B80D','#F16C00','#F22486']
class dfplt:
    def __init__(self,fignum,figsize):
        # Data Checks
        if fignum == 1:
            raise KeyError('Figure 1 is reserved for figure spacing graph. Please update fignum.')
        
        # Setup Default Parameters
        self.default_parameters()

        # Setup Lists for later use
        self.subplots = []
        self.axs = []
        self.figure_list = []


        # Setup Figure
        self.figsize = figsize
        self.fig = self.figure(fignum,figsize=figsize)

    def figure(self,fignum,figsize):
        # Figure size conversion: Points to Inches
        pt_to_inch = 1/72

        # Generate primary figure based on figure size
        fig = plt.figure(fignum,figsize=(pt_to_inch*figsize[0],pt_to_inch*figsize[1]))

        # Save figure number in setup
        self.figure_list.append(fig)

        return fig
    
    def load_csv(self,file_name):
        """
        Loads a CSV file and returns it as a pandas DataFrame.
        
        :param file_name: str, The path to the CSV file to load
        :return: DataFrame, Loaded pandas DataFrame
        """
        try:
            df = pd.read_csv(file_name)
            # print(f"CSV file '{file_name}' loaded successfully.")
            return df
        except FileNotFoundError:
            print(f"Error: The file '{file_name}' was not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

    def custom_cm(self,colors):
        cmap_name = 'custom_cmap'
        n_bins = 100
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

        return cm

    def default_parameters(self):
        # Text Parameters
        self.font = helvetica_font
        self.font_bold = helvetica_bold_font

        # Font size
        self.fontsize_label = 10
        self.fontsize_tick = 9
        self.fontsize_index = 12
    

        # Colors
        self.colors = colors
        self.color_bar = plt.cm.get_cmap('jet')
        # self.color_bar = ['#01B1CE','#3FB617','#F1B80D','#F16C00','#F22486']

    def plot(self, rect, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None,pad=0,projection=None):
        """
        Adds an axes to a given figure in a specified location with optional customization.
        
        Parameters:
            rect (list): A list specifying the axes position as [left, bottom, width, height].
            title (str): Title for the axes (optional).
            xlabel (str): Label for the x-axis (optional).
            ylabel (str): Label for the y-axis (optional).
            xlim (tuple): Limits for the x-axis (optional).
            ylim (tuple): Limits for the y-axis (optional).
            
        Returns:
            ax (matplotlib.axes.Axes): The created axes object.
        """

        # Add an axes to the figure at the specified location
        if projection=='polar':
            ax = self.fig.add_axes(rect,projection='polar')
        else:
            ax = self.fig.add_axes(rect)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', labelsize=10, pad=pad, length=6, direction='in')

        
        # Set title and labels if provided
        if title:
            ax.set_title(title, fontproperties=self.font_bold,fontsize=12,pad=pad)
        if xlabel:
            ax.set_xlabel(xlabel, fontproperties=self.font,fontsize=self.fontsize_label,labelpad=pad)
        if ylabel:
            ax.set_ylabel(ylabel, fontproperties=self.font,fontsize=self.fontsize_label,labelpad=pad)
        
        # Set axis limits if provided
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        # Remove the top and right boundaries
        if projection != 'polar':
            ax.spines['top'].set_visible(True)
            ax.spines['right'].set_visible(True)

        # Store the ax for reference
        self.axs.append(ax)

        # Return the axes to be used in a bigger plot later
        return ax
    
    def set_ticks(self,ax,xticks,yticks,length=None):
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        for label in ax.get_yticklabels():
            label.set_fontproperties(self.font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(self.font)
        # Set the font size for ticks
        if length != None:
            ax.tick_params(axis='x', labelsize=self.fontsize_tick,length =3)  # Change the font size
            ax.tick_params(axis='y', labelsize=self.fontsize_tick,length =3)  # Change the font size
        else:
            ax.tick_params(axis='x', labelsize=self.fontsize_tick)  # Change the font size
            ax.tick_params(axis='y', labelsize=self.fontsize_tick)  # Change the font size
        
        return
    
    def plot_png(self, rect, png_file, add_index=True):
        """
        Adds an axes to a given figure in a specified location and add a PNG file as an image.
        
        Parameters:
            rect (list): A list specifying the axes position as [left, bottom, width, height].
            png_file (str): Path to the PNG file.
            
        Returns:
            ax (matplotlib.axes.Axes): The created axes object.
        """
        if add_index:
            # Add the subplots to the class
            self.subplots.append(rect)

        # Add an axes to the figure at the specified location
        ax = self.fig.add_axes(rect)

        # Load the PNG image
        img = plt.imread(png_file)

        # Set the background of the axis to be transparent
        ax.patch.set_alpha(0.0)

        # Display the image using ax.imshow
        ax.imshow(img)

        # Remove all spines (the boundaries around the plot)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Remove the axes ticks
        ax.set_xticks([])
        ax.set_yticks([])

        # Store the ax for reference
        self.axs.append(ax)

        return ax
    
    def plot_pdf(self, rect, pdf_file,add_index=True):
        """
        Adds an axes to a given figure in a specified location and add an SVG file
        
        Parameters:
            rect (list): A list specifying the axes position as [left, bottom, width, height].
            
            
        Returns:
            ax (matplotlib.axes.Axes): The created axes object.
        """
        if add_index:
            # Add the subplots to the class
            self.subplots.append(rect)

        # Add an axes to the figure at the specified location
        ax = self.fig.add_axes(rect)

        # Convert PDF to image (first page)
        images = convert_from_path(pdf_file, dpi=1000)
        img = images[0]  # Get the first page

        # Display the image using ax.imshow
        ax.imshow(img)

        # Remove all spines (the boundaries around the plot)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Remove the axes ticks
        ax.set_xticks([])
        ax.set_yticks([])

        # Store the ax for reference
        self.axs.append(ax)

        return ax
    
    def add_cm_plot(self,ax,cm,norm,label='Strain (%)',labelpad=3,ticks=None):
        # Add a colorbar to show the relationship between line color and value
        # Create a ScalarMappable for the colorbar
        sm = ScalarMappable(cmap=cm, norm=norm)
        sm.set_array([])  # Dummy array for ScalarMappable

        # Add a colorbar
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(label, fontsize=self.fontsize_label, labelpad=labelpad,fontproperties=self.font)  # Set the font size of the label
        cbar.ax.tick_params(labelsize=self.fontsize_tick)  # Set the font size of the tick labels
        cbar.set_ticks(ticks)  # Set ticks every 5 units
    
        for label in cbar.ax.get_xticklabels():
                label.set_fontproperties(self.font)
        for label in cbar.ax.get_yticklabels():
            label.set_fontproperties(self.font)

        return cbar

    def set_tick_default(self,ax,pad=3):
        for label in ax.get_xticklabels():
            label.set_fontproperties(self.font)
        ax.tick_params(labelsize=self.fontsize_tick,pad=pad)  # Set the font size of the tick labels

    def plot_heatmap(self,data,ax,ticks,cb_label,pos,orientation='horizontal',NaN_colorbar=False,pad_cb=0,vmin=None,vmax=None,vertical=False):
        # Create a colormap
        custom_cmap = self.color_bar#plt.cm.colors.LinearSegmentedColormap.from_list('custom_cmap', self.color_bar, N=256)
        # colors = ['#1a1d5f', '#0064a1', '#00b5c3', '#48b612', '#eab800', '#f87e12', '#d90429']
        if NaN_colorbar:
            heat_map = custom_cmap
            # heat_map = plt.cm.colors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)
            newcolors = heat_map(np.linspace(0, 1, 256))
            white = np.array([1, 1, 1, 1])
            newcolors[0] = white  # Assuming the lower bound is 0
            newcolors[-1] = white  # Assuming the upper bound is 255
            custom_cmap = ListedColormap(newcolors)
        ax.axis('off')
        cax = ax.imshow(data, cmap=custom_cmap,aspect='auto',vmin=vmin,vmax=vmax)
        # Generate colorbar
        if vertical:
            cbar = plt.colorbar(cax, ax=ax,orientation='vertical',pad=pad_cb)
        else:
            cbar = plt.colorbar(cax, ax=ax,orientation='horizontal',pad=pad_cb)
        cbar.ax.set_position(pos) 
        # Set custom ticks
        cbar.set_ticks(ticks)
        
        # Set the label for the colorbar
        cbar.set_label(cb_label, fontsize=self.fontsize_label,fontproperties=self.font,labelpad=pad_cb)

        # Customize tick size and font
        cbar.ax.tick_params(labelsize=self.fontsize_tick)  # Set tick size
        cbar.ax.set_xticklabels(cbar.ax.get_xticks(), fontsize=self.fontsize_tick)  # Set tick font style

        for label in ax.get_xticklabels():
            label.set_fontproperties(self.font)
        # Set the font size for ticks
        ax.tick_params(axis='x', labelsize=self.fontsize_tick)  # Change the font size
        ax.tick_params(axis='y', labelsize=self.fontsize_tick)  # Change the font size
        
        
        return cbar
    
    def figure_grid(self):
        # Generate grid figure, always set as Figure 1
        self.fig_grid = self.figure(1,figsize=self.figsize)

        # if len(self.subplots) == 0:
        #     raise KeyError('No subplots have been saved. Unable to generate figure grid.')

        # Convert these relative values to absolute figure coordinates
        fig_width, fig_height = self.fig.get_size_inches()

        for rect in self.subplots:
            rect_width = fig_width * rect[2]
            rect_height = fig_height * rect[3]

            # Create the rectangle with absolute values
            rect_patch = patches.Rectangle((rect[0] * fig_width, rect[1] * fig_height), rect_width, rect_height,
                         linewidth=2, edgecolor='r', facecolor='none')
            
            # Add an axes to the figure at the specified location
            ax = self.fig_grid.add_axes(rect)
            # Remove ticks and labels
            ax.set_xticks([])  # Remove x ticks
            ax.set_yticks([])  # Remove y ticks
            ax.set_xticklabels([])  # Remove x tick labels
            ax.set_yticklabels([])  # Remove y tick labels

            # Add the rectangle to the plot
            ax.set_facecolor('lightblue')
            
            # Set the x and y axis limits based on figure size for proper scaling
            ax.set_xlim(0, fig_width)
            ax.set_ylim(0, fig_height)

            # Set equal scaling for better visualization
            # ax.set_aspect('equal')

        # Messaging
        print('Grid generated successfully. Saving file to "figure_grid.png".')
        plt.savefig('figure_grid.png')

    def add_index(self,convention='abc',padx=0,pady=0,font='bold'):
        if font=='bold':
            font = self.font_bold
        else:
            font=self.font

        if len(self.axs)==0:
            print('Note: Zero axes exist in plot, recheck code sequence.')
            return
        
        # Figure size conversion: Points to Inches
        pt_to_inch = 1/72

        # Marker labels
        markers = ['a', 'b', 'c', 'd','e','f','g']

        # Loop over the subplots
        for i, ax in enumerate(self.subplots):
            # Get the bounding box of the current axis in figure coordinates
            # bbox = ax.get_position()
            
            # Place the text in the top-left corner of each subplot
            # bbox.x0, bbox.y1 corresponds to the top-left of the axis
            self.fig.text(ax[0]+padx*pt_to_inch, ax[1]+ax[3]+pady*pt_to_inch, markers[i], fontsize=self.fontsize_index, fontweight='bold', ha='left', va='top')#,fontfamily=self.font_bold)
        
        # Return from function
        return

    def show(self,block=True):
        # Show the figures generated to the user
        plt.show(block=block)

        # Return after showing
        return
    
    def savefig(self,file_name,dpi,bbox_inches='tight'):

        # Save figure with tight layout
        plt.figure(2)
        self.fig.savefig(file_name,dpi=dpi)

            # Return after saving
        return
    
    # Create a function to crop the image into a smaller circle and draw the border
    def circular_label_image(self,image_name,ax, color,linewidth=2.25,rotate=False,padding_factor=0.05,file_type='.jpg'):
        cwd = os.getcwd()
        # Open the image
        img = Image.open(image_name+file_type)  # assuming PNG, change extension if needed
        size = img.size  # assuming image is square

        # Rotate the image by 90 degrees if rotate is True
        if rotate:
            img = img.rotate(90, expand=True)  # Rotate 90 degrees and expand canvas to fit

        # Calculate padding
        padding = int(size[0] * padding_factor)
        ellipse_box = (padding, padding, size[0] - padding, size[1] - padding)

        # Create a circular mask with smaller dimensions
        mask = Image.new('L', (size[0], size[1]), 0)
        draw = ImageDraw.Draw(mask)
        draw.ellipse(ellipse_box, fill=255)

        # Apply the circular mask
        circular_img = Image.new("RGBA", img.size)
        circular_img.paste(img, (0, 0), mask=mask)

        # Convert to numpy array for Matplotlib
        img_array = np.array(circular_img)

        # Draw the circular border using Matplotlib
        ax.imshow(img_array)

        # Create the circular border (smaller than the full image size)
        circle = plt.Circle((size[0] / 2, size[1] / 2), (size[0] / 2) - padding, color=color, fill=False, linewidth=linewidth)
        ax.add_patch(circle)

        # Hide axes and save the result
        ax.axis('off')
        return
    
def aspect_convert(locx,locy,figsize_sub,figsize):
    rect = [locx,locy,figsize_sub[0]/figsize[0],figsize_sub[1]/figsize[1]]
    return rect

if __name__ == "__main__":
    # Basic Testing of figure generator
    # Create a large main figure
    fig = dfplt(fignum=2,figsize=[240,240]) # Size in pt's

    # Create and add multiple axes to the main figure
    # The rect is in the form [left, bottom, width, height] (all values are relative to the figure size)

    ax1 = fig.plot([0.1, 0.7, 0.3, 0.2], xlabel="X", ylabel="Y", xlim=(0, 10), ylim=(0, 100))
    ax2 = fig.plot([0.5, 0.7, 0.4, 0.2], xlabel="X", ylabel="Y", xlim=(0, 10), ylim=(0, 50))
    ax3 = fig.plot([0.1, 0.4, 0.8, 0.2], xlabel="X", ylabel="Y", xlim=(0, 10), ylim=(0, 20))
    ax4 = fig.plot([0.3, 0.1, 0.5, 0.2], xlabel="X", ylabel="Y", xlim=(0, 10), ylim=(0, 5))

    # Plot sample data on each axes
    ax1.plot([1, 2, 3], [10, 40, 90], color='blue')
    ax2.plot([1, 2, 3], [5, 15, 45], color='green')
    ax3.plot([1, 2, 3], [2, 10, 18], color='red')
    ax4.plot([1, 2, 3], [1, 3, 4], color='purple')

    # Generate the grid spacing figure
    fig.figure_grid()

    # Show the main figure with all the axes
    fig.show()

    # Saving of figure
    fig.savefig('figure_1.pdf')