#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <tuple>
#include <utility>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xsort.hpp> // for argmin
#include <xtensor/xview.hpp> // for view
#include <xtensor/xbuilder.hpp> // for empty
#include <xtensor-blas/xlinalg.hpp> // for dot
#include <xtensor/xbuilder.hpp> // for xtuple
#include "_CFNode.h"
class _CFNode;

class Birch {
    // Implements the BIRCH clustering algorithm.

    // It is a memory-efficient, online-learning algorithm provided as an
    // alternative to :class:`MiniBatchKMeans`. It constructs a tree
    // data structure with the cluster centroids being read off the leaf.
    // These can be either the final cluster centroids or can be provided as input
    // to another clustering algorithm such as :class:`AgglomerativeClustering`.

    // Read more in the :ref:`User Guide <birch>`.

    // .. versionadded:: 0.16

    // Parameters
    // ----------
    // threshold : float, default=0.5
    //     The radius of the subcluster obtained by merging a new sample and the
    //     closest subcluster should be lesser than the threshold. Otherwise a new
    //     subcluster is started. Setting this value to be very low promotes
    //     splitting and vice-versa.

    // branching_factor : int, default=50
    //     Maximum number of CF subclusters in each node. If a new samples enters
    //     such that the number of subclusters exceed the branching_factor then
    //     that node is split into two nodes with the subclusters redistributed
    //     in each. The parent subcluster of that node is removed and two new
    //     subclusters are added as parents of the 2 split nodes.

    // n_clusters : int, instance of sklearn.cluster model or None, default=3
    //     Number of clusters after the final clustering step, which treats the
    //     subclusters from the leaves as new samples.

    //     - `None` : the final clustering step is not performed and the
    //       subclusters are returned as they are.

    //     - :mod:`sklearn.cluster` Estimator : If a model is provided, the model
    //       is fit treating the subclusters as new samples and the initial data
    //       is mapped to the label of the closest subcluster.

    //     - `int` : the model fit is :class:`AgglomerativeClustering` with
    //       `n_clusters` set to be equal to the int.

    // compute_labels : bool, default=True
    //     Whether or not to compute labels for each fit.

    // copy : bool, default=True
    //     Whether or not to make a copy of the given data. If set to False,
    //     the initial data will be overwritten.

    // Attributes
    // ----------
    // root_ : _CFNode
    //     Root of the CFTree.

    // dummy_leaf_ : _CFNode
    //     Start pointer to all the leaves.

    // subcluster_centers_ : ndarray
    //     Centroids of all subclusters read directly from the leaves.

    // subcluster_labels_ : ndarray
    //     Labels assigned to the centroids of the subclusters after
    //     they are clustered globally.

    // labels_ : ndarray of shape (n_samples,)
    //     Array of labels assigned to the input data.
    //     if partial_fit is used instead of fit, they are assigned to the
    //     last batch of data.

    // n_features_in_ : int
    //     Number of features seen during :term:`fit`.

    //     .. versionadded:: 0.24

    // feature_names_in_ : ndarray of shape (`n_features_in_`,)
    //     Names of features seen during :term:`fit`. Defined only when `X`
    //     has feature names that are all strings.

    //     .. versionadded:: 1.0

    // Notes
    // -----
    // The tree data structure consists of nodes with each node consisting of
    // a number of subclusters. The maximum number of subclusters in a node
    // is determined by the branching factor. Each subcluster maintains a
    // linear sum, mol_indices and the number of samples in that subcluster.
    // In addition, each subcluster can also have a node as its child, if the
    // subcluster is not a member of a leaf node.

    // For a new point entering the root, it is merged with the subcluster closest
    // to it and the linear sum, mol_indices and the number of samples of that
    // subcluster are updated. This is done recursively till the properties of
    // the leaf node are updated.

    // References
    // ----------
    // Original BIRCH
    // * Tian Zhang, Raghu Ramakrishnan, Maron Livny
    //   BIRCH: An efficient data clustering method for large databases.
    //   https://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf

    // * Roberto Perdisci
    //   JBirch - Java implementation of BIRCH clustering algorithm
    //   https://code.google.com/archive/p/jbirch

    public:
        float threshold;
        int branching_factor;
        int n_clusters;
        bool compute_labels;
        bool copy;
        int index_tracker;
        bool first_call;
        _CFNode* root_;
        _CFNode* dummy_leaf_;
        xt::xarray<float> subcluster_centers_;
        xt::xarray<float> subcluster_labels_;
        xt::xarray<float> labels_;
        int n_features_in_;
        xt::xarray<float> feature_names_in_;
        int _n_features_out;

        Birch(float threshold=0.5, int branching_factor=50, int n_clusters=3, bool compute_labels=true, bool copy=true);
        Birch fit(xt::xarray<float> X, void* y=NULL);
        Birch _fit(xt::xarray<float> X, bool partial);
        std::vector<_CFNode*> _get_leaves();
};