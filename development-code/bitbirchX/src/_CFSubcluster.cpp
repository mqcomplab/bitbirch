#include "_CFSubcluster.h"

_CFSubcluster::_CFSubcluster(xt::xarray<float> linear_sum, std::vector<int> mol_indices) {
    if (linear_sum.size() == 0) {
        this-> n_samples_ = 0;
        this-> centroid_ = this->linear_sum_ = 0;
        this->mol_indices = {};
    }
    else {
        this->n_samples_ = 1;
        this->centroid_ = this->linear_sum_ = linear_sum;
        this->mol_indices = mol_indices;
    }
    this->child_ = nullptr;
}

void _CFSubcluster::update(_CFSubcluster* subcluster) {
    this->n_samples_ += subcluster->n_samples_;
    this->linear_sum_ += subcluster->linear_sum_;
    this->mol_indices.insert(this->mol_indices.end(), subcluster->mol_indices.begin(), subcluster->mol_indices.end());
    this->centroid_ = calc_centroid(this->linear_sum_, this->n_samples_);
}

bool _CFSubcluster::merge_subcluster(_CFSubcluster* nominee_cluster, double threshold) {
    // Check if a cluster is worthy enough to be merged. If
    // yes then merge.
    xt::xarray<int> new_ls = this->linear_sum_ + nominee_cluster->linear_sum_;
    int new_n = this->n_samples_ + nominee_cluster->n_samples_;
    xt::xarray<float> new_centroid = calc_centroid(new_ls, new_n);
    std::vector<int> new_mol_indices = this->mol_indices;
    new_mol_indices.insert(new_mol_indices.end(), nominee_cluster->mol_indices.begin(), nominee_cluster->mol_indices.end());

    // TODO: Incorporate other criteria for merging

    float jt_radius = 1 - (jt_isim(new_ls + new_centroid, new_n + 1) * (new_n + 1) - jt_isim(new_ls, new_n) * (new_n - 1))/2; // Addition not to be taken literally

    if (jt_radius <= threshold) {
        this->n_samples_ = new_n;
        this->linear_sum_ = new_ls;
        this->centroid_ = new_centroid;
        this->mol_indices = new_mol_indices;
        return true;
    }
    return false;
}