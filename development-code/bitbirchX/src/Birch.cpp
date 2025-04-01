#include "Birch.h"

Birch::Birch(
    float threshold,
    int branching_factor,
    int n_clusters,
    bool compute_labels,
    bool copy
) {
    this->threshold = threshold;
    this->branching_factor = branching_factor;
    this->n_clusters = n_clusters;
    this->compute_labels = compute_labels;
    this->copy = copy;
    this->index_tracker = 0; // To support online updates
    this->first_call = true; // TODO: first_call wasn't properly updated between runs. Check this hack
}

Birch Birch::fit(const xt::xarray<float>& X, void* y) {
    // Build a CF Tree for the input data.

    // Parameters
    // ----------
    // X : {array-like, sparse matrix} of shape (n_samples, n_features)
    //     Input data.

    // y : Ignored
    //     Not used, present here for API consistency by convention.

    // Returns
    // -------
    // self
    //     Fitted estimator.

    // TODO: Add input verification

    return this->_fit(X, false);
}

Birch Birch::_fit(const xt::xarray<float>& X, bool partial) {
    threshold = this->threshold;
    branching_factor = this->branching_factor;

    int n_samples = X.shape(0);
    int n_features = X.shape(1);

    // If partial_fit is called for the first time or fit is called, we
    // start a new tree.
    if (this->first_call) {
        // The first root is the leaf. Manipulate this object throughout.
        this->root_ = new _CFNode(
            threshold,
            branching_factor,
            true,
            n_features
        );

        // To enable getting back subclustres.
        this->dummy_leaf_ = new _CFNode(
            threshold,
            branching_factor,
            true,
            n_features
        );
        this->dummy_leaf_->next_leaf_ = this->root_;
        this->root_->prev_leaf_ = this->dummy_leaf_;
    }

    // Cannot vectorize. Enough to convince to use cython.
    for (int i = 0; i < X.shape(0); i++) { 
        xt::xarray<int> sample = xt::row(X, i);
        _CFSubcluster* subcluster = new _CFSubcluster(sample, {this->index_tracker});
        bool split = this->root_->insert_cf_subcluster(subcluster);

        if (split) {
            std::pair<_CFSubcluster*, _CFSubcluster*> new_subclusters = _split_node(
                this->root_, threshold, branching_factor
            );
            delete this->root_;
            this->root_ = new _CFNode(
                threshold,
                branching_factor,
                false,
                n_features
            );
            this->root_->append_subcluster(new_subclusters.first);
            this->root_->append_subcluster(new_subclusters.second);
        }
        this->index_tracker += 1;
    }

    xt::xarray<float> centroids = this->_get_leaves()[0]->centroids_;
    bool skipFirst = true;
    for (_CFNode* leaf : this->_get_leaves()) {
        if (skipFirst) {
            skipFirst = false;
        }
        else {
            centroids = xt::concatenate(xt::xtuple(centroids, leaf->centroids_));
        }
    }
    this->subcluster_centers_ = centroids;
    this->_n_features_out = this->subcluster_centers_.shape(0);

    // TODO: Incorporate global_clustering option
    //this->_global_clustering(X);
    this->first_call = false;
    return *this;
}

std::vector<_CFNode*> Birch::_get_leaves() {
    //
    // Retrieve the leaves of the CF Node.
    //
    // Returns
    // -------
    // leaves : list of shape (n_leaves,)
    //  List of the leaf nodes.
    //
    _CFNode* leaf_ptr = this->dummy_leaf_->next_leaf_;
    std::vector<_CFNode*> leaves = {};
    while (leaf_ptr != nullptr) {
        leaves.push_back(leaf_ptr);
        leaf_ptr = leaf_ptr->next_leaf_;
    }
    return leaves;
}

Birch::~Birch() {
    // if (root_ != nullptr) {
    //     delete root_;
    //     root_ = nullptr;
    // }
    // if (dummy_leaf_ != nullptr) {
    //     delete dummy_leaf_;
    //     dummy_leaf_ = nullptr;
    // }
}