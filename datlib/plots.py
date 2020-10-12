#plots.py

# os module allows you to access commandline functions from python
import os
import pandas
# Math and data library
import numpy as np
import matplotlib.pyplot as plt

def plot_ts_scatter(df, s = 75, figsize = (40, 20), 
                    save_fig = False, pp = None):
   # Gather variables from df
    plot_vars = list(df.keys())
    for x in plot_vars:
        for y in plot_vars:
            if x != y:
                fig, ax = plt.subplots(figsize = figsize)
                # Create list of years from index
                # Years will be represented by color
                # Years will be the c value
                if "Year" not in df.keys():
                    # crate list from index
                    # convert each index vale to string
                    # only include first 4 characters, which is the y value
                    # Create an integer fro those characters
                    df["Year"] = [int(str(ind)[:4]) for ind in df.index] 
                df.plot.scatter(x = x, y = y, s = s, ax = ax, 
                                c = "Year", cmap = "viridis")
                
                # Turn the text on the x-axis so that it reads vertically
                ax.tick_params(axis='x', rotation=90)
                # Get rid of tick lines 
                ax.tick_params("both", length=0, which="both")
                # save image if PdfPages object was passed
                if save_fig:
                    try:
                        os.mkdir("plots")
                    except:
                        pass
                    # Identify directory to save figure
                    directory = "plots/" + x[:12] + " " + y[:12] + "c=Year" 
                    plt.savefig(directory + ".png")
                if pp != None: pp.savefig(fig, bbox_inches = "tight")
                    
                    
def plot_lines(df, linewidth = 1, figsize = (40,20), 
               legend = True, pp = None):
    fig, ax = plt.subplots(figsize = figsize)
    # If no secondary_y (axis), plot all variables at once
    df.plot.line(linewidth = linewidth, ax = ax, legend = legend)
    # Turn the text on the x-axis so that it reads vertically
    ax.tick_params(axis="x", rotation=90)
    # get rid of tick lines
    ax.tick_params("both", length=0, which = "both")
    
    vals = ax.get_yticks()
    vals = [int(x) for x in vals]
    ax.set_yticklabels(vals)
    
    # format image filename
    remove_chars = "[]:$'\\"
    filename = str(list(df.keys()))
    for char in remove_chars:
        filename = filename.replace(char, "")
    # avoid cutting off text
    plt.savefig(filename[:50] + "line.png",
               bbox_inches = "tight")
    if pp != None: pp.savefig(fig, box_inches = "tight")        


def plot_stacked_lines(df, plot_vars, linewidth = 1, figsize = (40,20), 
                       pp = None, total_var = False,
                       title = False):
    fig, ax = plt.subplots(figsize = figsize)
#    mpl_colors = ["C" + str(i) for i in range(11)]
    df[plot_vars].plot.area(stacked = True, linewidth = linewidth,
                            ax = ax)
    # change y vals from mil to tril
    
    if total_var != False:
        df[total_var].plot.line(linewidth = linewidth, ax = ax, c = "k",
              label = total_var, ls = "--")
    # place legend in top left corner of plot 
    # formal legend so that there are two columns of norms
    ax.legend(loc=2, ncol = 2)
    if title != False:
        plt.title(title)
         