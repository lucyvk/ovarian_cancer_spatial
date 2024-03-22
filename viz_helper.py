# Author: Lucy Van Kleunene
# Helper functions for spatial and network feature visualization

from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
import spatial_net_helper

def visualize_sample(fov_int,excluded,folder,points,types,show):
    
    type_colors = get_cell_type_colors()
    
    points_ex = np.array(points[fov_int])
    plot_colors = []
    for i in range(0,len(types[fov_int])):
        plot_colors.append(type_colors[types[fov_int][i]])

    fig = plt.figure()
    ax = fig.add_subplot()
    plt.title(f"Sample ID: {fov_int}")
    plt.figsize=(10,10)
    plt.xlim([0,1023])
    plt.ylim([0,1023])
    s = [1.2]*len(points_ex[:,0])
    plt.scatter(points_ex[:,0],points_ex[:,1],color=plot_colors,s=s)
    ax.set_aspect('equal', adjustable='box')
    plt.xticks([])
    plt.yticks([])

    if excluded:
        name = f'{folder}/dot_map_{fov_int}_excluded.png'
    else:
        name = f'{folder}/dot_map_{fov_int}.png'
    plt.savefig(name, dpi=1000, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
    
    # dimension assertions
    points_ex = np.array(points[fov_int])
    assert min(points_ex[:,0]) >= 0, "problem - min x dimension"
    assert max(points_ex[:,0]) < 1023, "problem - max x dimension"
    assert min(points_ex[:,1]) >= 0, "problem - min y dimension"
    assert max(points_ex[:,1]) < 1023, "problem - max y dimension"

    
def visualize_triangulation(fov_int,folder,points,types):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html
    # Unless you pass in the Qhull option “QJ”, Qhull does not guarantee that 
    # each input point appears as a vertex in the Delaunay triangulation. 
    # Omitted points are listed in the coplanar attribute.
    
    type_colors = get_cell_type_colors()
      
    points_ex = np.array(points[fov_int])
    plot_colors = []
    for i in range(0,len(types[fov_int])):
        plot_colors.append(type_colors[types[fov_int][i]])
    tri = Delaunay(points_ex, qhull_options="QJ")
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.title(f"Sample ID: {fov_int}")
    plt.figsize=(10,10)
    plt.xlim([0,1023])
    plt.ylim([0,1023])
    s = [1.2]*len(points_ex[:,0])
    plt.triplot(points_ex[:,0], points_ex[:,1], tri.simplices, linewidth=0.2, color='k')
    plt.scatter(points_ex[:,0],points_ex[:,1], color=plot_colors,s=s)
    ax.set_aspect('equal', adjustable='box')
    plt.xticks([])
    plt.yticks([])

    plt.savefig(f'{folder}/triangulation_{fov_int}.png', dpi=1000, bbox_inches='tight')
    plt.show()

    
def visualize_spatial_network(fov_int,folder,points,types,th):
    
    type_colors = get_cell_type_colors()
    
    spatial_net, pos_dict, removed, _ = spatial_net_helper.create_spatial_network(fov_int,points,types,th,False)
    plot_colors = []
    for i in range(0,len(types[fov_int])):
        plot_colors.append(type_colors[types[fov_int][i]])
        
    fig = plt.figure()
    plt.figsize=(10,10)
    ax = fig.add_subplot()

    num_cells = len(types[fov_int])

    # re-visualize from the networkx object
    nx.draw_networkx(spatial_net,pos=pos_dict,node_color=plot_colors, node_size=1.2, width=0.2, with_labels=False, ax=ax)
    plt.xlim([0,1023])
    plt.ylim([0,1023])
    ax.set_aspect('equal', adjustable='box')
    plt.xticks([])
    plt.yticks([])
    plt.title(f"Sample ID: {fov_int}") #, {num_cells} cells, distance threshold: {th}, {removed} edges removed")
    
    plt.savefig(f'{folder}/triangulation_network_trimmed_{th}_id_{fov_int}.png',dpi=1000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

def visualize_spatial_network_regions(fov_int,folder,points,types,th):
    
    type_colors = get_cell_type_colors()
    
    spatial_net, pos_dict, removed,_ = spatial_net_helper.create_spatial_network(fov_int,points,types,th,True)
    plot_colors = []
    for i in range(0,len(types[fov_int])):
        plot_colors.append(type_colors[types[fov_int][i]])
    
    fig = plt.figure()
    plt.figsize=(10,10)
    ax = fig.add_subplot()
    num_cells = len(types[fov_int])
    # re-visualize from the networkx object
    nx.draw_networkx(spatial_net,pos=pos_dict,node_color=plot_colors, node_size=1.2, width=0.2, with_labels=False,ax=ax)
    plt.xlim([0,1023])
    plt.ylim([0,1023])
    ax.set_aspect('equal', adjustable='box')
    plt.xticks([])
    plt.yticks([])
    plt.title(f"Sample ID: {fov_int}")#, {num_cells} cells, distance threshold: {th}")
    
    plt.savefig(f'{folder}/network_regions_{th}_id_{fov_int}.png', dpi=1000,bbox_inches='tight')    
    plt.show()


def get_cell_type_colors():
    cell_type_colors = {'CD8+ T cells': '#6a3d9a', 'CD4+ T cells': '#b15928', 'NK/NKT': '#ff7f00', 'CD56+CD45-': '#1f78b4',
                        'Vascular endothelial cells': '#ebe652', 'Lymphatic endothelial cells': '#fb9a99', 'Fibroblast': '#33a02c',
                        'Tumor': '#e31a1c', 'Unidentified': '#a6cee3', 'B cells': '#bc80bd', 'M2 macrophages': '#cab2d6',
                        'M1 macrophages': '#b2df8a', 'Dendritic cells': '#d9d9d9', 'CD163+ cells': '#eb52d4', 'Monocytes': '#000000',
                        'CD11c_low immune': '#4d4f4e', 'Other immune': '#fdbf6f', 'Neuroepithelial cells': '#2e3785',
                        'Non-leukocyte derived neural cells': '#15f4ee', 'CD11b+ epithelial': '#bf00ff', 'Neutrophils': '#aaf001',
                        'CD11c+ epithelial': '#fe347e', 'CD11b_low Neutrophils': '#ff9966', 'HLADR+': '#00ff00'}
    return cell_type_colors
    