/* Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
          Joel Nothman <joel.nothman@gmail.com>
 License: BSD 3 clause */

#include "abstract_birch.h"

// xt::xarray<float> find_centroid(xt::xarray<float> X) {
//     // find centroid of matrix X
//     // return centroid;
// }

// float isim(xt::xarray<float> X) {
//     // calculate iSIM for matrix X
//     // return isim; 
// } 

// xt::xarray<float> sim_search(xt::xarray<float> point, xt::xarray<float> X) { 
//     // Find similarities between a point and all the rows of matrix X
//     // return similarities; 
// } 

float jt_isim(xt::xarray<float> c_total, int n_objects) {
    // iSIM Tanimoto calculation
    float sum_kq = xt::sum(c_total)();
    float sum_kqsq = xt::linalg::dot(c_total, c_total)();
    float a = (sum_kqsq - sum_kq)/2;
    return a/(a + n_objects * sum_kq - sum_kqsq);
}

std::tuple<std::pair<int, int>, xt::xarray<float>, xt::xarray<float>> max_separation(xt::xarray<float> X) {
    // Finds two objects in X that are very separated
    // This is an approximation (is not guaranteed to find
    // the two absolutely most separated objects), but it is
    // a very robust O(N) approximation
    
    // Algorithm:
    // a) Find centroid of X
    // b) mol1 is the molecule most distant from the centroid
    // c) mol2 is the molecule most distant from mol1
    
    // Returns
    // -------
    // (mol1, mol2) : (int, int)
    //                indices of mol1 and mol2
    // 1 - sims_mol1 : xt::xarray<int>
    //                Distances to mol1
    // 1 - sims_mol2: xt::xarray<int>
    //                Distances to mol2
    // These are needed for node1_dist and node2_dist in _split_node

    int n_samples = X.shape(0);
    xt::xarray<float> linear_sum = xt::sum(X, 0);
    xt::xarray<float> med = calc_centroid(linear_sum, n_samples);
    xt::xarray<float> pop_counts = xt::sum(X, 1); // the number of ones for every row
    xt::xarray<float> a_med = xt::linalg::dot(X, med); // numerator of the tanimoto between the centroid and every data point (1D vector)
    xt::xarray<float> sims_med = a_med / (xt::sum(med) + pop_counts - a_med);
    int mol1 = xt::argmin(sims_med)[0];
    xt::xarray<float> a_mol1 = xt::linalg::dot(X, xt::view(X, mol1, xt::all()));
    xt::xarray<float> sims_mol1 = a_mol1 / (xt::view(pop_counts, mol1, xt::all()) + pop_counts - a_mol1);
    int mol2 = xt::argmin(sims_mol1)[0];
    xt::xarray<float> a_mol2 = xt::linalg::dot(X, xt::view(X, mol2, xt::all()));
    xt::xarray<float> sims_mol2 = a_mol2 / (xt::view(pop_counts, mol2, xt::all()) + pop_counts - a_mol2);
    return make_tuple(std::make_pair(mol1, mol2), 1 - sims_mol1, 1 - sims_mol2);

    // // abstract_birch implementation
    // xt::xarray<float> centroid = find_centroid(X);
    // xt::xarray<int> mol1 = xt::argmin(sim_search(centroid, X));
    // xt::xarray<int> mol2 = xt::argmin(sim_search(mol1, X));
    // return make_tuple(std::make_pair(mol1, mol2), 1 - sim_search(mol1, X), 1 - sim_search(mol2, X));
}

xt::xarray<float> calc_centroid(xt::xarray<float> linear_sum, int n_samples, std::string sim_index) { // linear_sum is the same as c_total
    // Calculates centroid
    if (sim_index == "JT") {
        return xt::floor(linear_sum / n_samples + 0.5);
    }
}

// xt::xarray<float> _iterate_sparse_X(xt::xarray<float> X) {
//     // This little hack returns a densified row when iterating over a sparse
//     // matrix, instead of constructing a sparse  matrix for every row that is
//     // expensive.
//     //
//     int n_samples = X.shape(0);
//     xt::xarray<int> X_indices = X.indices; //these are from scipy.sparse import csr_matrix
//     xt::xarray<float> X_Data = X.data();
//     xt::xarray<int> X_indptr = X.indptr;

//     xt::xarray<float> result;
//     for (int i = 0; i < n_samples; i++) {
//         auto row = xt::zeros<float>(X.shape(1));
//         int startptr = X_indptr[i];
//         int endptr = X_indptr[i + 1];
//         xt::xarray<int> nonzero_indices = X_indices(X_indices.begin() + startptr, X_indices.begin() + endptr);
//         xt::xarray<float> row[nonzero_indices] = X_data[X_indices.begin() + startptr, X_indices.begin() + endptr];
//         result.push_back(row);
//     }
//     return result;
// }

std::pair<_CFSubcluster*, _CFSubcluster*> _split_node(_CFNode* node, double threshold, int branching_factor) {
    // The node has to be split if there is no place for a new subcluster
    // in the node.
    // 1. Two empty nodes and two empty subclusters are initialized.
    // 2. The pair of distant subclusters are found.
    // 3. The properties of the empty subclusters and nodes are updated
    //    according to the nearest distance between the subclusters to the
    //    pair of distant subclusters.
    // 4. The two nodes are set as children to the two subclusters.
    _CFSubcluster* new_subcluster1 = new _CFSubcluster();
    _CFSubcluster* new_subcluster2 = new _CFSubcluster();
    _CFNode* new_node1 = new _CFNode(
        threshold=threshold, 
        branching_factor=branching_factor, 
        node->is_leaf, 
        node->n_features //----maybe private these attributes and use getters?
        // dtype=node->init_centroids_.dtype //----dtype is from numpy
    );
    _CFNode* new_node2 = new _CFNode(
        threshold=threshold, 
        branching_factor=branching_factor, 
        node->is_leaf, 
        node->n_features //----maybe private these attributes and use getters?
        // dtype=node->init_centroids_.dtype //----dtype is from numpy
    );
    new_subcluster1->child_ = new_node1;
    new_subcluster2->child_ = new_node2;

    if (node->is_leaf) {
        if (node->prev_leaf_ != nullptr) {
            node->prev_leaf_->next_leaf_ = new_node1; 
        }
        new_node1->prev_leaf_ = node->prev_leaf_;
        new_node1->next_leaf_ = new_node2;
        new_node2->prev_leaf_ = new_node1;
        new_node2->next_leaf_ = node->next_leaf_;
        if (node->next_leaf_ != nullptr) {
            node->next_leaf_->prev_leaf_ = new_node2;
        }
    }

    std::tuple<std::pair<xt::xarray<int>, xt::xarray<int>>, xt::xarray<float>, xt::xarray<float>> result = max_separation(node->centroids_);
    std::pair<xt::xarray<int>, xt::xarray<int>> farthest_idx = std::get<0>(result);
    xt::xarray<float> node1_dist = std::get<1>(result);
    xt::xarray<float> node2_dist = std::get<2>(result);

    xt::xarray<float> node1_closer = node1_dist < node2_dist;
    // make sure node1 is closest to itself even if all distances are equal.
    // This can only happen when all node.centroids_ are duplicates leading to all
    // distances between centroids being zero.
    node1_closer[farthest_idx.first] = true;

    for (int i = 0; i < node->subclusters_.size(); i++) { 
        _CFSubcluster* subcluster = node->subclusters_[i];
        if (node1_closer[i]) {
            new_node1->append_subcluster(subcluster);
            new_subcluster1->_CFSubcluster::update(subcluster);
        }
        else {
            new_node2->append_subcluster(subcluster);
            new_subcluster2->update(subcluster);
        }
    }
    return std::make_pair(new_subcluster1, new_subcluster2);
}

// float cluster_radius(xt::xarray<float> cluster) { 
//     // Calculate cluster radius
//     //return radius
// } 