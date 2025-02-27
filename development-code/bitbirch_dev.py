# BitBIRCH is an open-source clustering module based on iSIM
#
# THIS IS A DEVELOPMENT VERSION OF THE CODE
# THIS INCLUDES SOME CHANGES NOT PRESENT IN THE STABLE RELEASE
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
#
# BitBIRCH License: LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.en.html#license-text
#
# Please, cite the BitBIRCH paper: https://www.biorxiv.org/content/10.1101/2024.08.10.607459v1
#
### Part of the tree-management code was derived from https://scikit-learn.org/stable/modules/generated/sklearn.cluster.Birch.html
### Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
###          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
###          Joel Nothman <joel.nothman@gmail.com>
### License: BSD 3 clause

import time
import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin

def jt_distances(X):
    """Calculates the matrix of Tanimoto distances
    This is VERY inefficient
    """
    dist = []
    for i in X:
        dist.append([])
        si = np.sum(i)
        for j in X:
            a = np.dot(i,j)
            dist[-1].append(1 - a/(si + np.sum(j) - a))
    dist = np.array(dist)
    return dist


def jt_isim(c_total, n_objects):
    """iSIM Tanimoto calculation"""
    sum_kq = np.sum(c_total)
    sum_kqsq = np.dot(c_total, c_total)
    a = (sum_kqsq - sum_kq)/2

    return a/(a + n_objects * sum_kq - sum_kqsq)

def max_separation(X):
    """Finds two objects in X that are very separated
    This is an approximation (is not guaranteed to find
    the two absolutely most separated objects), but it is
    a very robust O(N) approximation
    
    Algorithm:
    a) Find centroid of X
    b) mol1 is the molecule most distant from the centroid
    c) mol2 is the molecule most distant from mol1
    
    Returns
    -------
    (mol1, mol2) : (int, int)
                   indices of mol1 and mol2
    1 - sims_mol1 : np.ndarray
                   Distances to mol1
    1 - sims_mol2: np.ndarray
                   Distances to mol2
    These are needed for node1_dist and node2_dist in _split_node
    """
    n_samples = len(X)
    linear_sum = np.sum(X, axis = 0)
    med = calc_centroid(linear_sum, n_samples)
    pop_counts = np.sum(X, axis = 1)

    a_med = np.dot(X, med)
    sims_med = a_med / (np.sum(med) + pop_counts - a_med)

    mol1 = np.argmin(sims_med)
    a_mol1 = np.dot(X, X[mol1])
    sims_mol1 = a_mol1 / (pop_counts[mol1] + pop_counts - a_mol1)
    
    mol2 = np.argmin(sims_mol1)
    a_mol2 = np.dot(X, X[mol2])
    sims_mol2 = a_mol2 / (pop_counts[mol2] + pop_counts - a_mol2)
    
    return (mol1, mol2), sims_mol1, sims_mol2


def calc_centroid(linear_sum, n_samples, sim_index = 'JT'):
    """Calculates centroid"""
    return np.where(linear_sum >= n_samples * 0.5, 1, 0)

def _iterate_sparse_X(X):
    """This little hack returns a densified row when iterating over a sparse
    matrix, instead of constructing a sparse matrix for every row that is
    expensive.
    """
    n_samples, n_features = X.shape
    X_indices = X.indices
    X_data = X.data
    X_indptr = X.indptr

    for i in range(n_samples):
        row = np.zeros(n_features)
        startptr, endptr = X_indptr[i], X_indptr[i + 1]
        nonzero_indices = X_indices[startptr:endptr]
        row[nonzero_indices] = X_data[startptr:endptr]
        yield row


def _split_node(node, threshold, branching_factor):
    """The node has to be split if there is no place for a new subcluster
    in the node.
    1. Two empty nodes and two empty subclusters are initialized.
    2. The pair of distant subclusters are found.
    3. The properties of the empty subclusters and nodes are updated
       according to the nearest distance between the subclusters to the
       pair of distant subclusters.
    4. The two nodes are set as children to the two subclusters.
    """
    new_subcluster1 = _BFSubcluster()
    new_subcluster2 = _BFSubcluster()
    new_node1 = _BFNode(
        threshold=threshold,
        branching_factor=branching_factor,
        is_leaf=node.is_leaf,
        n_features=node.n_features,
        dtype=node.init_centroids_.dtype,
    )
    new_node2 = _BFNode(
        threshold=threshold,
        branching_factor=branching_factor,
        is_leaf=node.is_leaf,
        n_features=node.n_features,
        dtype=node.init_centroids_.dtype,
    )
    new_subcluster1.child_ = new_node1
    new_subcluster2.child_ = new_node2

    if node.is_leaf:
        if node.prev_leaf_ is not None:
            node.prev_leaf_.next_leaf_ = new_node1
        new_node1.prev_leaf_ = node.prev_leaf_
        new_node1.next_leaf_ = new_node2
        new_node2.prev_leaf_ = new_node1
        new_node2.next_leaf_ = node.next_leaf_
        if node.next_leaf_ is not None:
            node.next_leaf_.prev_leaf_ = new_node2
    
    
    # O(N) implementation of max separation
    farthest_idx, node1_dist, node2_dist = max_separation(node.centroids_)    
    
    node1_closer = node1_dist > node2_dist
    node1_closer[farthest_idx[0]] = True

    for idx, subcluster in enumerate(node.subclusters_):
        if node1_closer[idx]:
            new_node1.append_subcluster(subcluster)
            new_subcluster1.update(subcluster)
        else:
            new_node2.append_subcluster(subcluster)
            new_subcluster2.update(subcluster)
        
    new_subcluster1.centroid_ = calc_centroid(new_subcluster1.linear_sum_, new_subcluster1.n_samples_)
    new_subcluster2.centroid_ = calc_centroid(new_subcluster2.linear_sum_, new_subcluster2.n_samples_)

    return new_subcluster1, new_subcluster2


class _BFNode:
    """Each node in a BFTree is called a BFNode.

    The BFNode can have a maximum of branching_factor
    number of BFSubclusters.

    Parameters
    ----------
    threshold : float
        Threshold needed for a new subcluster to enter a BFSubcluster.

    branching_factor : int
        Maximum number of BF subclusters in each node.

    is_leaf : bool
        We need to know if the BFNode is a leaf or not, in order to
        retrieve the final subclusters.

    n_features : int
        The number of features.

    Attributes
    ----------
    subclusters_ : list
        List of subclusters for a particular BFNode.

    prev_leaf_ : _BFNode
        Useful only if is_leaf is True.

    next_leaf_ : _BFNode
        next_leaf. Useful only if is_leaf is True.
        the final subclusters.

    init_centroids_ : ndarray of shape (branching_factor + 1, n_features)
        Manipulate ``init_centroids_`` throughout rather than centroids_ since
        the centroids are just a view of the ``init_centroids_`` .

    centroids_ : ndarray of shape (branching_factor + 1, n_features)
        View of ``init_centroids_``.

    """

    def __init__(self, *, threshold, branching_factor, is_leaf, n_features, dtype):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.is_leaf = is_leaf
        self.n_features = n_features

        # The list of subclusters, centroids and squared norms
        # to manipulate throughout.
        self.subclusters_ = []
        self.init_centroids_ = np.zeros((branching_factor + 1, n_features), dtype=dtype)
        self.prev_leaf_ = None
        self.next_leaf_ = None

        self.centroids_sum=np.zeros((1, branching_factor + 1),dtype="int64")

    def append_subcluster(self, subcluster):
        n_samples = len(self.subclusters_)

        self.subclusters_.append(subcluster)
        self.init_centroids_[n_samples] = subcluster.centroid_
        
        # Keep centroids as views. In this way
        # if we change init_centroids, it is sufficient
        self.centroids_ = self.init_centroids_[: n_samples + 1, :]

        self.centroids_sum[0, n_samples] = np.sum(subcluster.centroid_)

        
    def update_split_subclusters(self, subcluster_ind, new_subcluster1, new_subcluster2):
        """Remove a subcluster from a node and update it with the
        split subclusters.
        """
        ind = subcluster_ind
        self.subclusters_[ind] = new_subcluster1
        self.init_centroids_[ind] = new_subcluster1.centroid_
        self.centroids_[ind] = new_subcluster1.centroid_
        self.append_subcluster(new_subcluster2)

        self.centroids_sum[0, ind] = np.sum(new_subcluster1.centroid_)

    def insert_bf_subcluster(self, subcluster):
        """Insert a new subcluster into the node."""
        if not self.subclusters_:
            self.append_subcluster(subcluster)
            return False

        threshold = self.threshold
        branching_factor = self.branching_factor
        # We need to find the closest subcluster among all the
        # subclusters so that we can insert our new subcluster.

        a = np.dot(self.centroids_, subcluster.centroid_)
        sim_matrix = self.centroids_sum[0][:len(self.centroids_)] / a # inverse 
        closest_index = np.argmin(sim_matrix) # inverse
        closest_subcluster = self.subclusters_[closest_index]

        # If the subcluster has a child, we need a recursive strategy.
        if closest_subcluster.child_ is not None:
            split_child = closest_subcluster.child_.insert_bf_subcluster(subcluster)

            if not split_child:
                # If it is determined that the child need not be split, we
                # can just update the closest_subcluster
                closest_subcluster.update(subcluster)
                temp=self.subclusters_[closest_index].centroid_
                self.init_centroids_[closest_index] = temp
                self.centroids_[closest_index] = temp
                self.centroids_sum[0, closest_index] = np.sum(temp)
                return False

            # things not too good. we need to redistribute the subclusters in
            # our child node, and add a new subcluster in the parent
            # subcluster to accommodate the new child.
            else:
                new_subcluster1, new_subcluster2 = _split_node(
                    closest_subcluster.child_,
                    threshold,
                    branching_factor
                )
                self.update_split_subclusters(
                    closest_index, new_subcluster1, new_subcluster2
                )

                if len(self.subclusters_) > self.branching_factor:
                    return True
                return False

        # good to go!
        else:
            merged = closest_subcluster.merge_subcluster(subcluster, self.threshold)
            if merged:
                temp=closest_subcluster.centroid_
 
                self.centroids_[closest_index]=temp             
                self.init_centroids_[closest_index] = temp
                self.centroids_sum[0, closest_index] = np.sum(temp)

                return False

            # not close to any other subclusters, and we still
            # have space, so add.
            elif len(self.subclusters_) < self.branching_factor:
                self.append_subcluster(subcluster)
                return False

            # We do not have enough space nor is it closer to an
            # other subcluster. We need to split.
            else:
                self.append_subcluster(subcluster)
                return True


class _BFSubcluster:
    """Each subcluster in a BFNode is called a BFSubcluster.

    A BFSubcluster can have a BFNode has its child.

    Parameters
    ----------
    linear_sum : ndarray of shape (n_features,), default=None
        Sample. This is kept optional to allow initialization of empty
        subclusters.

    Attributes
    ----------
    n_samples_ : int
        Number of samples that belong to each subcluster.

    linear_sum_ : ndarray
        Linear sum of all the samples in a subcluster. Prevents holding
        all sample data in memory.

    centroid_ : ndarray of shape (branching_factor + 1, n_features)
        Centroid of the subcluster. Prevent recomputing of centroids when
        ``BFNode.centroids_`` is called.
    
    mol_indices : list, default=[]
        List of indices of molecules included in the subclustergiven cluster.

    child_ : _BFNode
        Child Node of the subcluster. Once a given _BFNode is set as the child
        of the _BFNode, it is set to ``self.child_``.
    """

    def __init__(self, *, linear_sum = None, mol_indices = []):
        if linear_sum is None:
            self.n_samples_ = 0
            self.centroid_ = self.linear_sum_ = 0
            self.mol_indices = []
        else:
            self.n_samples_ = 1
            self.centroid_ = self.linear_sum_ = linear_sum
            self.mol_indices = mol_indices
        self.child_ = None

    def update(self, subcluster):
        self.n_samples_ += subcluster.n_samples_
        self.linear_sum_ += subcluster.linear_sum_
        self.mol_indices += subcluster.mol_indices

    def merge_subcluster(self, nominee_cluster, threshold):
        """Check if a cluster is worthy enough to be merged. If
        yes then merge.
        """
        new_ls = self.linear_sum_ + nominee_cluster.linear_sum_
        new_n = self.n_samples_ + nominee_cluster.n_samples_
        new_centroid = calc_centroid(new_ls, new_n)
        
        # TODO: Incorporate other criteria for merging
        # multiply 2 instead of div

        ls_cent=new_ls+new_centroid
        jt_radius = jt_isim(ls_cent, new_n + 1) * (new_n + 1) - jt_isim(new_ls, new_n) * (new_n - 1)
        
        if jt_radius >= threshold*2:
            (
                self.n_samples_,
                self.linear_sum_,
                self.centroid_,
                self.mol_indices,
            ) = (new_n, new_ls, new_centroid, self.mol_indices + nominee_cluster.mol_indices)

            return True
        return False


class BitBirch():
    """Implements the BitBIRCH clustering algorithm.
    
    BitBIRCH paper: 

    Memory- and time-efficient, online-learning algorithm.
    It constructs a tree data structure with the cluster centroids being read off the leaf.
    
    Parameters
    ----------
    threshold : float, default=0.5
        The similarity radius of the subcluster obtained by merging a new sample and the
        closest subcluster should be greater than the threshold. Otherwise a new
        subcluster is started. Setting this value to be very low promotes
        splitting and vice-versa.

    branching_factor : int, default=50
        Maximum number of BF subclusters in each node. If a new samples enters
        such that the number of subclusters exceed the branching_factor then
        that node is split into two nodes with the subclusters redistributed
        in each. The parent subcluster of that node is removed and two new
        subclusters are added as parents of the 2 split nodes.

    Attributes
    ----------
    root_ : _BFNode
        Root of the BFTree.

    dummy_leaf_ : _BFNode
        Start pointer to all the leaves.

    subcluster_centers_ : ndarray
        Centroids of all subclusters read directly from the leaves.

    subcluster_labels_ : ndarray
        Labels assigned to the centroids of the subclusters after
        they are clustered globally.
    
    labels_ : ndarray of shape (n_samples,)
        Array of labels assigned to the input data.

    Notes
    -----
    The tree data structure consists of nodes with each node consisting of
    a number of subclusters. The maximum number of subclusters in a node
    is determined by the branching factor. Each subcluster maintains a
    linear sum, mol_indices and the number of samples in that subcluster.
    In addition, each subcluster can also have a node as its child, if the
    subcluster is not a member of a leaf node.

    For a new point entering the root, it is merged with the subcluster closest
    to it and the linear sum, mol_indices and the number of samples of that
    subcluster are updated. This is done recursively till the properties of
    the leaf node are updated.
    """


    def __init__(
        self,
        *,
        threshold=0.5,
        branching_factor=50,
        n_clusters=3,
        compute_labels=True,
        copy=True,
        perform_clustering=False,
    ):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.n_clusters = n_clusters
        self.compute_labels = compute_labels
        self.copy = copy
        self.index_tracker = 0
        self.first_call = True
        self.perform_clustering = perform_clustering

    def fit(self, X, y=None):
        """
        Build a BF Tree for the input data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input data.

        y : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        self
            Fitted estimator.
        """

        # TODO: Add input verification
        return self._fit(X, partial=False)

    def _fit(self, X, partial):
        threshold = self.threshold
        branching_factor = self.branching_factor

        n_samples, n_features = X.shape
        

        # If partial_fit is called for the first time or fit is called, we
        # start a new tree.
        if self.first_call:
            # The first root is the leaf. Manipulate this object throughout.
            self.root_ = _BFNode(
                threshold=threshold,
                branching_factor=branching_factor,
                is_leaf=True,
                n_features=n_features,
                dtype=X.dtype,
            )
    
            # To enable getting back subclusters.
            self.dummy_leaf_ = _BFNode(
                threshold=threshold,
                branching_factor=branching_factor,
                is_leaf=True,
                n_features=n_features,
                dtype=X.dtype,
            )
            self.dummy_leaf_.next_leaf_ = self.root_
            self.root_.prev_leaf_ = self.dummy_leaf_

        # Cannot vectorize. Enough to convince to use cython.
        if not sparse.issparse(X):
            iter_func = iter
        else:
            iter_func = _iterate_sparse_X

        for sample in iter_func(X):

            subcluster = _BFSubcluster(linear_sum=sample, mol_indices = [self.index_tracker])

            split = self.root_.insert_bf_subcluster(subcluster)

            if split:
                new_subcluster1, new_subcluster2 = _split_node(
                    self.root_, threshold, branching_factor
                )
                del self.root_
                self.root_ = _BFNode(
                    threshold=threshold,
                    branching_factor=branching_factor,
                    is_leaf=False,
                    n_features=n_features,
                    dtype=X.dtype,
                    
                )
                self.root_.append_subcluster(new_subcluster1)
                self.root_.append_subcluster(new_subcluster2)
            self.index_tracker += 1

        centroids = np.concatenate([leaf.centroids_ for leaf in self._get_leaves()])
        self.subcluster_centers_ = centroids
        self._n_features_out = self.subcluster_centers_.shape[0]
        
        if(self.perform_clustering):
            self._global_clustering(X)
        self.first_call = False
        return self

    def _get_leaves(self):
        """
        Retrieve the leaves of the BF Node.

        Returns
        -------
        leaves : list of shape (n_leaves,)
            List of the leaf nodes.
        """
        leaf_ptr = self.dummy_leaf_.next_leaf_
        leaves = []
        while leaf_ptr is not None:
            leaves.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf_
        return leaves
    
    def retrieveVal(self):
        print()

    def _global_clustering(self, X):
        clusters = self.n_clusters
        centroids = self.subcluster_centers_
        compute_labels = (X is not None) and self.compute_labels

        if isinstance(clusters, int):
            km = KMeans(n_clusters=clusters)
            self.subcluster_labels_ = km.fit_predict(centroids)
        else:
            # argument is None (skip global clustering)
            self.subcluster_labels_ = np.arange(len(centroids))
            return

        if compute_labels:
            argmin = pairwise_distances_argmin(X, centroids)
            self.labels_ = self.subcluster_labels_[argmin]