# BitBIRCH is an open-source clustering module based on iSIM
#
# Please, cite the BitBIRCH paper: https://www.biorxiv.org/content/10.1101/2024.08.10.607459v1
#
# BitBIRCH is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# BitBIRCH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# BitBIRCH authors (PYTHON): Ramon Alain Miranda Quintana <ramirandaq@gmail.com>, <quintana@chem.ufl.edu>
#                            Vicky (Vic) Jung <jungvicky@ufl.edu>
#                            Kenneth Lopez Perez <klopezperez@chem.ufl.edu>
#                            Kate Huddleston <kdavis2@chem.ufl.edu>
#
# BitBIRCH License: LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.en.html#license-text

import numpy as np
import copy

# when deleting from init_centroids, we need to make sure that the elements are shifted
# to the left to fill the space and make sure the empty space on the right is 0 because
# init_centroids is an array with a pre-defined size
def shift_delete(array, index):
    for i in range(index, len(array)-1):
        array[i]=array[i+1]
    array[-1]=0
    
def recursivelyTraverse(node, children):
    if node is None:
        return
    
    for i in node.subclusters_:
        children.append(i)
        recursivelyTraverse(i.child_, children)

def pSubPath(subcluster, pSubs, brc):
    # figure out the path of subclusters from leaf level to root

    if subcluster is None:
        return 
    
    if subcluster.parent_ is None:
        pSubs.append(brc.root_.subclusters_.index(subcluster))
    else:
        pSubs.append(subcluster.parent_.child_.subclusters_.index(subcluster))
    
    pSubPath(subcluster.parent_, pSubs, brc)

# in order to preserve the insertion order of pNodes, make sure python is updated to 3.7+
def pNodePath(pSubs, pNodes, root):
    temp=root

    for i in pSubs[:0:-1]:
        index=temp.subclusters_.index(i)
        pNodes[temp]=index
        temp=temp.subclusters_[index].child_


def calc_centroid(linear_sum, n_samples):
    """Calculates centroid
    
    Parameters
    ----------
    
    linear_sum : np.ndarray
                 Sum of the elements column-wise
    n_samples : int
                Number of samples
                
    Returns
    -------
    centroid : np.ndarray
               Centroid fingerprints of the given set
    """
    return np.where(linear_sum >= n_samples * 0.5 , 1, 0)

def lazyPrune(leaf_index, subc_index, brc):
    #from bitbirch_tools import bitbirch_double_link_tolerance as bitbirch
    # path from leaf level to root, indices are marked 
    pSubs=[]
    pSubPath(brc._get_leaves()[leaf_index].subclusters_[subc_index], pSubs, brc)

    # traverse all the way down to the leaf node
    node=brc.root_
    for index in pSubs[:0:-1]:
        node=node.subclusters_[index].child_

    
    # information needed to update all preceding subclusters
    cluster_ls = node.subclusters_[subc_index].linear_sum_
    cluster_n = node.subclusters_[subc_index].n_samples_
    cluster_mol_indices = node.subclusters_[subc_index].mol_indices
    
    '''# IMPORTANT: tree must have at least three levels, at least that is what is assumed when inputting dataset
    parent_sub = node.subclusters_[0].parent_       # needed to traverse node upwards
    parent_node = parent_sub.parent_.child_         # needed to reassign leaf nodes'''

    parent_sub = node.subclusters_[0].parent_
    if parent_sub.parent_ is None:
        parent_node = brc.root_
    else:
        parent_node = parent_sub.parent_.child_
    
    # eliminate the biggest subcluster from the leaf node
    node.subclusters_.pop(subc_index)
    node.centroids_ = np.delete(node.centroids_, subc_index, 0)
    shift_delete(node.init_centroids_, subc_index)

    if len(node.subclusters_) == 0:
        parent_sub.child_ = None

    # rearrange leaf pointers if need to

    # if the leaf node disappears, we assign parent node to replace leaf node
    num_childs=0
    for i in parent_node.subclusters_:
        if i.child_ is not None:
            num_childs+=1
    if num_childs == 0:
        parent_node.is_leaf = True
        if node.prev_leaf_ is not None:
            node.prev_leaf_.next_leaf_ = parent_node
        parent_node.prev_leaf_ = node.prev_leaf_
        parent_node.next_leaf_ = node.next_leaf_
        if node.next_leaf_ is not None:
            node.next_leaf_.prev_leaf_ = parent_node

    # leaf node disappears and there are no reassignments
    elif len(node.subclusters_) == 0:
        if node.prev_leaf_ is not None:
            node.prev_leaf_.next_leaf_ = node.next_leaf_
        if node.next_leaf_ is not None:
            node.next_leaf_.prev_leaf_ = node.prev_leaf_


    if parent_sub.parent_ is None:
        node = brc.root_
    else:
        node = parent_sub.parent_.child_ 

  
    # here, we update all the preceding subclusters all the way up to the root
    for index in pSubs[1:]:
        sub=node.subclusters_[index]
 
        sub.linear_sum_ -= cluster_ls
        sub.n_samples_ -= cluster_n
        sub.centroid_ = calc_centroid(sub.linear_sum_, sub.n_samples_)
        sub.mol_indices = [x for x in sub.mol_indices if x not in cluster_mol_indices]
        
        node.centroids_[index] = sub.centroid_
        node.init_centroids_[index] = sub.centroid_

        if sub.parent_ is None:
            break
        elif sub.parent_.parent_ is None:
            node=brc.root_
        else:
            node=sub.parent_.parent_.child_

    return cluster_mol_indices

