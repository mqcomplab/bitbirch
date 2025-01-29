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
#include <xtensor-blas/xlinalg.hpp> //for dot
#include "_CFSubcluster.h"
#include "abstract_birch.h"
class _CFSubcluster;

class _CFNode {
    // Each node in a CFTree is called a CFNode.

    // The CFNode can have a maximum of branching_factor
    // number of CFSubclusters.

    // Parameters
    // ----------
    // threshold : float
    //     Threshold needed for a new subcluster to enter a CFSubcluster.

    // branching_factor : int
    //     Maximum number of CF subclusters in each node.

    // is_leaf : bool
    //     We need to know if the CFNode is a leaf or not, in order to
    //     retrieve the final subclusters.

    // n_features : int
    //     The number of features.

    // Attributes
    // ----------
    // subclusters_ : list
    //     List of subclusters for a particular CFNode.

    // prev_leaf_ : _CFNode
    //     Useful only if is_leaf is True.

    // next_leaf_ : _CFNode
    //     next_leaf. Useful only if is_leaf is True.
    //     the final subclusters.

    // init_centroids_ : ndarray of shape (branching_factor + 1, n_features)
    //     Manipulate ``init_centroids_`` throughout rather than centroids_ since
    //     the centroids are just a view of the ``init_centroids_`` .

    // centroids_ : ndarray of shape (branching_factor + 1, n_features)
    //     View of ``init_centroids_``.

    public:
        float threshold;
        int branching_factor;
        bool is_leaf;
        int n_features;
        std::vector<_CFSubcluster*> subclusters_;
        _CFNode* prev_leaf_;
        _CFNode* next_leaf_;
        xt::xarray<int> init_centroids_;
        xt::xarray<int> centroids_;
        int n_samples;

        _CFNode(double threshold, int branching_factor, bool is_leaf, int n_features);
        void append_subcluster(_CFSubcluster* subcluster);
        void update_split_subclusters(_CFSubcluster* subcluster, _CFSubcluster* new_subcluster1, _CFSubcluster* new_subcluster2);
        bool insert_cf_subcluster(_CFSubcluster* subcluster);
}; 