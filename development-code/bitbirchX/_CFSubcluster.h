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
#include <xtensor/xsort.hpp> //for argmin
#include <xtensor/xview.hpp> // for view
#include <xtensor/xbuilder.hpp> // for empty
#include <xtensor-blas/xlinalg.hpp> //for dot
#include "_CFNode.h"
class _CFNode;

class _CFSubcluster {
    // Each subcluster in a CFNode is called a CFSubcluster.

    // A CFSubcluster can have a CFNode has its child.

    // Parameters
    // ----------
    // linear_sum : ndarray of shape (n_features,), default=None
    //     Sample. This is kept optional to allow initialization of empty
    //     subclusters.

    // Attributes
    // ----------
    // n_samples_ : int
    //     Number of samples that belong to each subcluster.

    // linear_sum_ : ndarray
    //     Linear sum of all the samples in a subcluster. Prevents holding
    //     all sample data in memory.

    // centroid_ : ndarray of shape (branching_factor + 1, n_features)
    //     Centroid of the subcluster. Prevent recomputing of centroids when
    //     ``CFNode.centroids_`` is called.
    
    // mol_indices : list, default=[]
    //     List of indices of molecules included in the subclustergiven cluster.

    // child_ : _CFNode
    //     Child Node of the subcluster. Once a given _CFNode is set as the child
    //     of the _CFNode, it is set to ``self.child_``.

    public:
        xt::xarray<float> linear_sum;
        int n_samples_;
        xt::xarray<float> linear_sum_;
        xt::xarray<float> centroid_;
        std::vector<int> mol_indices;
        _CFNode* child_;

        _CFSubcluster(xt::xarray<float> linear_sum = {}, std::vector<int> mol_indices = {});
        void update(_CFSubcluster* subcluster);
        bool merge_subcluster(_CFSubcluster* nominee_cluster, double threshold);
};