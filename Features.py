import numpy as np
from neuron_morphology.swc_io import morphology_from_swc

def get_distance_from_axon_origin_to_first_branch(morph):
    """
    Finds euclidean distance between axon origination and first axon branch node
    Inputs:
    morph: a neuron_morphology Morphology object
    
    Returns:
    distance: float
              euclidean distance between axon origination and first axon branch node
              
    first_axon_branch: Morphology node, dict
                       First branching node in axon tree
                       
    """
    first_axon_node = morph.get_node_by_types([2])[0]
    my_node = morph.get_node_by_types([2])[0]
    while len(morph.get_children(my_node)) < 2:
        if morph.get_children(my_node) == []:
            #There is not much axon detected in this cell, only one unbranching segment 
            return np.nan, my_node
        else:
            my_node = morph.get_children(my_node)[0]
        
    first_axon_branch = my_node
    
    distance = euclidean_dist_from_nodes(first_axon_node,first_axon_branch)
    return distance,first_axon_branch

def path_distance_from_nodes(ancestor,child,morph):
    return morph.get_segment_length([ancestor,child])
        
    
def euclidean_dist_from_nodes(n1,n2):
    return (((n1['x']-n2['x'])**2)  + ((n1['y']-n2['y'])**2)  + ((n1['z']-n2['z'])**2))**0.5
    

def get_path_distance_between_branch_nodes(root_branch_node,order,morph,output):
    """
    Function to find path distance between branch nodes.
    Root branch node = top branch node, 2nd value return from get_distance_from_axon_origin_to_first_branch
    order = how many branch nodes down the tree to explore
    morph = morphology object
    output = []. Initialize as an empty list and gets passed along through the recursive calls
    """
    if order >0:
        empty_list = []
        next_branch_nodes = find_next_branch_recursion(root_branch_node,morph,empty_list)
        for nod in next_branch_nodes:
            if morph.get_children(nod) != []:
                path_dist_btwn_branches = path_distance_from_nodes(root_branch_node,nod,morph)
                output.append(path_dist_btwn_branches)
                get_path_distance_between_branch_nodes(nod,order-1,morph,output)

    return output
    
def find_next_branch_recursion(tree_node,morph,ap_list):   
    """
    Recursive function to find the next set of branching nodes. 
    Inputs:
    tree_node: a morphology node who has more than one child
    morph: morphology object
    ap_list starts as an empty list 
    
    Note this does not claim that the next returned nodes will have children that branch
    """
    next_branch_node = []
    if morph.get_children(tree_node) != []:
        for my_node in morph.get_children(tree_node):
            if len(morph.get_children(my_node)) < 2:
                find_next_branch_recursion(my_node,morph,ap_list)
            else:
                ap_list.append(my_node)
    return ap_list


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    num = np.dot(v1,v2)
    denom = np.linalg.norm(v1)*np.linalg.norm(v2)

    return np.arccos(num/denom) 


def get_theta_between_branch_point_vectors(root_branch_node,order,morph,output):  
    """
    Recursive function to find angle between branching nodes
    Handles trifurcations as well by taking average of 3 angles
    
    Root_branch_node = top branch node to begin analysis. 2nd return value from get_distance_from_axon_origin_to_first_branch
    order = how far down the tree to go
    """
    if order >0:
        empty_list = []
        next_branch_nodes = find_next_branch_recursion(root_branch_node,morph,empty_list)
        root_coords = np.array([root_branch_node['x'],root_branch_node['y'],root_branch_node['z']])
        #Create vectors from root point to branch points
        parent_child_vector_dict = {}
        for no in next_branch_nodes:
            parent_child_vector_dict[no['id']] = root_coords-np.array([no['x'],no['y'],no['z']])
            
        parent_child_vector_dict_array = list(parent_child_vector_dict.values())
         
        if len(next_branch_nodes) > 2:
            theta1 = angle_between(parent_child_vector_dict_array[0],parent_child_vector_dict_array[1]) #angle between those lines
            theta2 = angle_between(parent_child_vector_dict_array[0],parent_child_vector_dict_array[2]) #angle between those lines
            theta3 = angle_between(parent_child_vector_dict_array[1],parent_child_vector_dict_array[2]) #angle between those lines
            theta = np.nanmean(np.asarray([theta1,theta2,theta3]))
            output.append(theta)
            print('trifurcation')
            
        elif len(next_branch_nodes) == 1:
            get_theta_between_branch_point_vectors(next_branch_nodes[0],order-1,morph,output)
            
        elif len(next_branch_nodes) == 0:
            bad_node = next_branch_nodes
               
        else:

            theta = angle_between(parent_child_vector_dict_array[0],parent_child_vector_dict_array[1]) #angle between those lines
            output.append(theta)          
            for nod in next_branch_nodes:
                if morph.get_children(nod) != []:
                    get_theta_between_branch_point_vectors(nod,order-1,morph,output)
        
    return output
