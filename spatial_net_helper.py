# Author: Lucy Van Kleunene
# Helper functions for spatial and network feature generation

import numpy as np
from scipy.spatial import Delaunay
import networkx as nx
from scipy.spatial import distance
from scipy.spatial import KDTree
import viz_helper

# Median (spatial) distance from {type_to} cells to cells of each other cell type in a sample (fov_int)
def generate_median_distances(fov_int, points, types, type_to):
    points_ex = points[fov_int]
    
    # The type_to is missing from the sample 
    if type_to not in types[fov_int]:
        return None
    
    # Construct KD Trees for nearest neighbor look ups for all points
    type_pts = {} # cell type : [list of (x,y) coordinates for cells of this type]
    for k in range(0,len(points_ex)):
        curr_type = types[fov_int][k]
        if curr_type not in type_pts:
            type_pts[curr_type] = [points_ex[k]]
        else:
            type_pts[curr_type].append(points_ex[k])
    kd_trees = {} # cell type: KD Tree 
    for curr_type in type_pts.keys():
        if curr_type != type_to:
            kd_trees[curr_type] = KDTree(type_pts[curr_type])
            
    nearest_distances = {}
    # Loop through type_to cells
    for cell in type_pts[type_to]:
        # For each other cell type
        for cell_type in kd_trees.keys():
            # Find the nearest neighbor of this cell type (KD tree makes this faster)
            nearest = kd_trees[cell_type].query(cell,k=1) 
            if cell_type in nearest_distances:
                nearest_distances[cell_type].append(nearest[0])
            else:
                nearest_distances[cell_type] = [nearest[0]]
    
    for cell_type in nearest_distances.keys():
        nearest_distances[cell_type] = np.median(nearest_distances[cell_type])
    return(nearest_distances)

def create_spatial_network(fov_int,points,types,th,regions):
    
    unique_types = set(viz_helper.get_cell_type_colors().keys())
    
    pos_dict = {}
    
    points_ex = np.array(points[fov_int]) # 2D points just for this sample    
    tri = Delaunay(points_ex, qhull_options="QJ") # triangulation
    
    attr_dict = {}
    
    # NetworkX Graph object
    spatial_net = nx.Graph()
    
    # loop through verticies and make dictionary mapping node id to position + add to network w/ class attribute
    for k in range(0,len(points_ex[:,0])):
        pos_dict[k] = list(points_ex[k,:]) # node id mapped to 2D position vector
        types_dict = {} # node id mapped to binary dictionary for cell type
        for key in unique_types: # For every unique cell type
            types_dict[key] = 0 # the attribute value for this node for that cell type is 0
        types_dict[types[fov_int][k]] = 1 # except for its true cell type
        types_dict["cell type"] = types[fov_int][k]
        attr_dict[k] = types_dict
        spatial_net.add_node(k)
        
    nx.set_node_attributes(spatial_net,attr_dict) # attributes are the binary cell type dictionaries
    local_distance_dict = {} # save all the distances

    removed = 0 
    indptr, indices = tri.vertex_neighbor_vertices
    # loop through verticies
    for k in range(0,len(points_ex[:,0])):
        # get neighbors and add edges 
        neighbors = indices[indptr[k]:indptr[k+1]]
        for n in neighbors:
            #check if the distance is less than a certain threshold and exclude
            cell_type1 = types[fov_int][k]
            cell_type2 = types[fov_int][n] # neighbor cell type 
            curr_distance = distance.euclidean(pos_dict[k],pos_dict[n]) # distance to spatial neighbor
            if curr_distance < th: # added if under trimming threshold
                if (regions and cell_type1 == cell_type2) or not regions: # only connect cells of the same type, optionally
                    spatial_net.add_edge(k,n)
                if (cell_type1, cell_type2) in local_distance_dict:
                    local_distance_dict[(cell_type1, cell_type2)].append(curr_distance)
                elif (cell_type2, cell_type1) in local_distance_dict: # check reversed version as well 
                    local_distance_dict[(cell_type2, cell_type1)].append(curr_distance)
                else:
                    local_distance_dict[(cell_type1, cell_type2)] = [curr_distance]
            else:
                removed += 1
    
    return spatial_net, pos_dict, removed, local_distance_dict

def contact_enrichment_score(fov_int, net, type_to, points_ex, types_ex):
    
    # The type_to is missing from the sample 
    if type_to not in types_ex:
        return None
    
    focal_cells = [] # list of ids of cells of type_to (focal cell type)
    label_list = [] # ordered list of cell labels for non-focal cells
    uni_types = [] 
    for k in range(0,len(points_ex)):
        curr_type = types_ex[k] # get type
        if curr_type == type_to:
            focal_cells.append(k)
        else:
            label_list.append(curr_type)
    
    uni_types = set(types_ex)
    uni_types.remove(type_to) # unique cell types not including focal cell type 
    label_list = np.array(label_list)
    
    # Look at real - For the non-shuffled version:
    # Loop through focal cells and look at neighbors
    real_freqs = {} # cell type : number of contacts with focal cell type 
    for cell_type in uni_types:
        real_freqs[cell_type] = 0 
    for focal_cell in focal_cells:
        neighbors = net.neighbors(focal_cell)
        for n in neighbors:
            cell_type2 = types_ex[n] 
            if cell_type2 != type_to:
                real_freqs[cell_type2] = real_freqs[cell_type2] + 1
        
    # Generate a null distribution of cell type contact frequencies 
    null_freqs = {} # {cell type: [list of 1000 frequencies]}
    for i in range(0,1000): 
        # Shuffle all of the cell labels except for the tumor cells
        np.random.shuffle(label_list)
        
        # Assign the new labels
        new_labels = {}
        x = 0 
        for k in range(0,len(points_ex)):
            if types_ex[k] != type_to:
                # Re-assign non-tumor cell types a label from the shuffled list 
                new_labels[k] = label_list[x]
                x+=1
            else:
                new_labels[k] = type_to
                
        assert x == len(label_list), " list length doesn't match "
            
        # Loop through tumor cells and look at neighbors
        temp_freqs = {}
        for cell_type in uni_types:
            temp_freqs[cell_type] = 0 
        for focal_cell in focal_cells:
            neighbors = net.neighbors(focal_cell)
            for n in neighbors:
                cell_type2 = new_labels[n]
                if cell_type2 != type_to:
                    temp_freqs[cell_type2] = temp_freqs[cell_type2] + 1
        
        for cell_type in temp_freqs:
            if cell_type not in null_freqs:
                null_freqs[cell_type] = [temp_freqs[cell_type]]
            else:
                null_freqs[cell_type].append(temp_freqs[cell_type])
        
    scores = {}
    for type_key in null_freqs:
        
        # Calculate Z score comparing real score to the null distribution
        curr_mean = np.mean(null_freqs[type_key])
        curr_std = np.std(null_freqs[type_key])
        if curr_std != 0:
            z_score = (real_freqs[type_key] - curr_mean)/curr_std
        else:
            z_score = 0
        scores[type_key] = z_score
    
    return(scores)
