import numpy as np
import matplotlib.pyplot as plt
def cartes_to_polar (x,y):
    r = np.sqrt(x**2+y**2)
    if x!=0:
        theta = np.arctan(y/x)
    else:
        theta = np.pi/2
    return r,theta
def draw_filtration_circles(elements, labels=None,outline = True,cmap = 'Blues'):
    """
    Draw co-centered, filled circles representing the number of elements at each step of the filtration.
    
    Args:
    elements (list of int): A list of the number of elements at each filtration step. 
                            The first element is the initial number of elements.
    labels (list of str, optional): A list of labels for each step. The number of elements will 
                                    automatically be appended to each label. Defaults to 'Step i' if not provided.
    
    The function draws circles with radii proportional to the number of elements left after each step,
    and fills them with different colors from a colormap.
    """
    if labels!=None:
        assert len(elements)==len(labels), "lengths and labels must be the same length"
    # Ensure that the input list is in decreasing order
    elements = sorted(elements, reverse=True)
    
    # Normalize the element sizes to be proportional to the circle radii
    max_elements = elements[0]
    radii = [elem / max_elements for elem in elements]
    
    # Get a colormap (viridis is an eye-pleasing colormap)
    cmap = plt.get_cmap(cmap)
    
    # Set up the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    # Remove spines and labels
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Draw each circle with a different fill color from the colormap
    for i, radius in enumerate(radii):
        color = cmap(i / len(radii))  # Get color from the colormap
        circle = plt.Circle((0, 0), radius, color=color, fill=True)
        ax.add_artist(circle)
    
    if outline:
        circle = plt.Circle((0, 0), radii[0], color='black', fill=False)
        ax.add_artist(circle)
        
    # Set limits to ensure all circles are visible
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    
    # Set up legend labels, append number of elements to each label
    if labels is None:
        labels = [f'Step {i+1}: {int(elements[i])} elements' for i in range(len(elements))]
    else:
        labels = [f'{labels[i]}: {int(elements[i])} elements' for i in range(len(labels))]

    # Add a legend outside of the plot to avoid covering the circles
    handles = [plt.Line2D([0], [0], color=cmap(i / len(radii)), lw=4) for i in range(len(radii))]
    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    return fig,ax

from matplotlib.patches import Wedge

def draw_filtration_quarter_circles(elements, labels=None, outline = True):
    """
    Draw quarter (first quadrant) circles representing the number of elements at each step of the filtration.
    
    Args:
    elements (list of int): A list of the number of elements at each filtration step. 
                            The first element is the initial number of elements.
    labels (list of str, optional): A list of labels for each step. The number of elements will 
                                    automatically be appended to each label. Defaults to 'Step i' if not provided.
    
    The function draws quarter circles with radii proportional to the number of elements left after each step,
    and fills them with different colors from the 'Blues' colormap.
    """
    # Ensure that the input list is in decreasing order
    elements = sorted(elements, reverse=True)
    
    # Normalize the element sizes to be proportional to the circle radii
    max_elements = elements[0]
    radii = [elem / max_elements for elem in elements]
    
    # Use the 'Blues' colormap
    cmap = plt.get_cmap('Blues')
    
    # Set up the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    # Remove x and y axis labels
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Keep the bottom and left spines, but remove top and right
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    
    # Draw each quarter circle (90 degree wedges) with a different fill color from the colormap
    for i, radius in enumerate(radii):
        color = cmap(i / len(radii))  # Get color from the colormap
        # Wedge starts at 0 degrees and ends at 90 degrees (first quadrant)
        quarter_circle = Wedge((0, 0), radius, 0, 90, color=color, fill=True)
        ax.add_artist(quarter_circle)
    if outline:
        quarter_circle = Wedge((0, 0), radii[0], 0, 90, color='black', fill=False)
        ax.add_artist(quarter_circle)

    # Set limits to ensure all quarter circles are visible
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    
    # Set up legend labels, append number of elements to each label
    if labels is None:
        labels = [f'Step {i+1}: {elements[i]} elements' for i in range(len(elements))]
    else:
        labels = [f'{labels[i]}: {elements[i]} elements' for i in range(len(labels))]

    # Add a legend outside of the plot to avoid covering the quarter circles
    handles = [plt.Line2D([0], [0], color=cmap(i / len(radii)), lw=4) for i in range(len(radii))]
    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    return fig, ax

import sys
import os

from Combo_CC_Base import Combo_CC
from Combo_plot_base import Combo_plot
import sympy as sp
def draw_combos_phase_portraits (cms,par_dicts,ncols = 4,verbose = False,
                                 ax_size_width = 5, ax_size_height = 5,
                                 x_lims = [-0.2,2],y_lims = [-0.2,2]):
    nrows = int(len(cms)/ncols)+int(len(cms)%ncols>0)
    fig,axs = plt.subplots(nrows = nrows, ncols = ncols,figsize = (int(ax_size_width*ncols),int(ax_size_height*nrows)))
    for i,combo in enumerate(cms):
        if verbose:
            print(f"plotting combo {i+1}/{len(cms)}")
        if nrows>1:
            ax = axs[int(i/ncols),i%ncols]
        elif ncols>1:
            ax = axs[i]
        else:
            ax = axs
        
        combo.find_nullclines()
        if type(par_dicts)==list:
            cm_pl = Combo_plot(combo,par_dicts[i])
        else:
            cm_pl = Combo_plot(combo,par_dicts)
        if type(x_lims[0])==list and type(y_lims[0])==list:
            cm_pl.plot_phase_portrait(0,30,30,x_lims[i],y_lims[i],colors = ['tab:blue','tab:orange'],fig = fig, ax = ax)
        elif type(x_lims[0])==list and type(y_lims[0])!=list:
            cm_pl.plot_phase_portrait(0,30,30,x_lims[i],y_lims,colors = ['tab:blue','tab:orange'],fig = fig, ax = ax)
        elif type(x_lims[0])!=list and type(y_lims[0])==list:
            cm_pl.plot_phase_portrait(0,30,30,x_lims,y_lims[i],colors = ['tab:blue','tab:orange'],fig = fig, ax = ax)
        elif type(x_lims[0])!=list and type(y_lims[0])!=list:
            cm_pl.plot_phase_portrait(0,30,30,x_lims,y_lims,colors = ['tab:blue','tab:orange'],fig = fig, ax = ax)
        
        p_latex = sp.latex(combo.P)
        q_latex = sp.latex(combo.Q)

        textstr = "$ \\dot{x}"+f" = {p_latex}$ \n" +"$\\dot{y} ="+f" {q_latex}$"
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.95, 0.05, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='bottom', horizontalalignment='right', bbox=props)
    fig.tight_layout()
    return fig,axs
import itertools
import pandas as pd

def generate_combinations(lists, k_values,list_names = None):
    # Ensure that the number of lists matches the number of k values
    assert len(lists) == len(k_values), "The number of lists must match the number of k_values."

    # Generate all combinations for each list
    combinations = []
    for lst, k in zip(lists, k_values):
        combinations.append(list(itertools.combinations(lst, k)))

    # Generate the Cartesian product of the combinations from each list
    all_combinations = list(itertools.product(*combinations))
    
    # Flatten the combinations (since each combination from a list is a tuple)
    flattened_combinations = [tuple(itertools.chain.from_iterable(comb)) for comb in all_combinations]

    # Create column names dynamically
    columns = []
    for i, k in enumerate(k_values):
        if list_names == None:
            columns.extend([f'List_{i+1}_item_{j+1}' for j in range(k)])
        else:
            columns.extend([f'{list_names[i]}_item_{j+1}' for j in range(k)])

    # Convert the result to a pandas DataFrame
    df = pd.DataFrame(flattened_combinations, columns=columns)

    return df

import matplotlib.patches as patches
def draw_circuit_from_cm(cm:Combo_CC,fig = None, ax = None,xlim = [-0.25,3.25],ylim =[0,2],
                         markersize = 20,d = 0.15,
                         input_node = False, input_to = None):
    if fig == None and ax==None:
        fig, ax = plt.subplots(figsize=(6, 6))
    positions = {"X": (0,1),"Y": (3,1)} 
    if input_node:
        positions ['u'] = (1.5,1.75)
    for node, (x, y) in positions.items():
        ax.plot(x, y, 'o', markersize=markersize, color='lightblue',zorder = 1)
        ax.text(x, y, node, ha='center', va='center', fontsize=10, fontweight='bold')
    def draw_arrow (start,end,effect,type='str',linestyle = '-'):
        arrowstyle = '->' if effect == "+" else '-['  # solid for positive, dashed for negative
        if type=='str':
            connecstyle = 'arc3'
        elif type=='circ':
            connecstyle = f'arc3, rad={.5*0.15/d}'
        elif type=='circ_ccw':
            connecstyle = f'arc3, rad={-.5*0.15/d}'
        elif type=='full_circ':
            connecstyle = f'arc3, rad={3*0.15/d}'
        elif type=='full_circ_ccw':
            connecstyle = f'arc3, rad={-3*0.15/d}'
        kw = dict(arrowstyle=arrowstyle, color="k",linestyle=linestyle,)
        a3 = patches.FancyArrowPatch(start, end, mutation_scale=10,
                             connectionstyle=connecstyle, **kw,zorder = 0)
        ax.add_patch(a3)
    for term in cm.prod_x_terms:
        if term=='1':
            draw_arrow((0,2-d),(0,1+d),'+')
        if term=='y':
            draw_arrow((0,2-d),(0,1+d),'+')
            draw_arrow((3-d,1),(d,1.5),'+')
        if term=='x':
            draw_arrow((d,1-d),(d,1+d),'+','full_circ')
    for term in cm.rem_x_terms:
        if term=='1':
            draw_arrow((0,1-d),(0,d),'+')
        if term=='x':
            draw_arrow((0,1-d),(0,d),'+')
            draw_arrow((-d,1),(-d,0.5),'+','circ')
        if term=='y':
            draw_arrow((3-d,1),(d,0.5),'+')
        if term=='x^2':
            draw_arrow((d,1-d),(d,1+d),'+','full_circ')
            draw_arrow((0,1+d),(0.75,1+d),'-','circ_ccw')
        if term=='x*y':
            draw_arrow((d,1-d),(d,1+d),'+','full_circ')
            draw_arrow((3-d,1),(1.3,1),'-')
        if term=='y^2':
            draw_arrow((0,2-d),(0,1+d),'+')
            draw_arrow((3-d,1),(d,0.5),'+',linestyle='--')
    for term in cm.prod_y_terms:
        if term=='1':
            draw_arrow((3,2-d),(3,1+d),'+')
        if term=='y':
            draw_arrow((3-d,1-d),(3-d,1+d),'+','full_circ_ccw')
        if term=='x':
            draw_arrow((3,2-d),(3,1+d),'+')
            draw_arrow((d,1),(3-d,1.5),'+')
    for term in cm.rem_y_terms:
        if term==1:
            draw_arrow((3,1-d),(3,d),'+')
        if term=='y':
            draw_arrow((3,1-d),(3,d),'+')
            draw_arrow((3+d,1),(3+d,0.5),'+','circ_ccw')
        if term=='x':
            draw_arrow((3,1-d),(3,d),'+')
            draw_arrow((d,1),(3-d,0.5),'+')
        if term=='y^2':
            draw_arrow((3-d,1-d),(3-d,1+d),'+','full_circ_ccw')
            draw_arrow((3,1+d),(2.25,1+d),'-','circ')
        if term=='x*y':
            draw_arrow((3-d,1-d),(3-d,1+d),'+','full_circ_ccw')
            draw_arrow((d,1),(2-2*d,1),'-')
        if term=='x^2':
            draw_arrow((3,1-d),(3,d),'+')
            draw_arrow((d,1),(3-d,0.5),'+',linestyle = '--')
    if input_node:
        if input_to=='x':
            draw_arrow((1.5-d,1.75-d),(0,1+d),'+')
        elif input_to=='y':
            draw_arrow((1.5+d,1.75-d),(3,1+d),'+')


    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_axis_off()
    return fig, ax

