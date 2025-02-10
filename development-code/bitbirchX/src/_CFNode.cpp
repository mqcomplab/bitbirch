#include "_CFNode.h"

_CFNode::_CFNode(double threshold, int branching_factor, bool is_leaf, int n_features) {
    this->threshold = threshold;
    this->branching_factor = branching_factor;
    this->is_leaf = is_leaf;
    this->n_features = n_features;

    // The list of subclusters, centroids and squared norms
    // to manipulate throughout.
    this->subclusters_ = {};
    this->init_centroids_ = xt::zeros<float>({branching_factor + 1, n_features}); 
    this->prev_leaf_ = nullptr;
    this->next_leaf_ = nullptr;
}

void _CFNode::append_subcluster(_CFSubcluster* subcluster) {
    this->n_samples = this->subclusters_.size();
    this->subclusters_.push_back(subcluster);
    xt::row(this->init_centroids_, n_samples) = subcluster->centroid_;

    // Keep centroids as views. In this way
    // if we change init_centroids, it is sufficient
    this->centroids_ = xt::view(this->init_centroids_, xt::range(0, n_samples + 1), xt::all());
}

void _CFNode::update_split_subclusters(_CFSubcluster* subcluster, _CFSubcluster* new_subcluster1, _CFSubcluster* new_subcluster2) {
    // Remove a subcluster from a node and update it with the 
    // split subclusters.
    int ind = find(this->subclusters_.begin(), this->subclusters_.end(), subcluster) - this->subclusters_.begin();
    this->subclusters_[ind] = new_subcluster1;
    xt::row(this->init_centroids_, ind) = new_subcluster1->centroid_;
    this->centroids_ = xt::view(this->init_centroids_, xt::range(0, this->n_samples + 1), xt::all());
    this->append_subcluster(new_subcluster2);
}

bool _CFNode::insert_cf_subcluster(_CFSubcluster* subcluster) {
    //Insert a new subcluster into the node.
    if (this->subclusters_.empty()) {
        this->append_subcluster(subcluster);
        return false;
    }

    threshold = this->threshold;
    branching_factor = this->branching_factor;
    // We need to find the closest subcluster among all the
    // subclusters so that we can insert our new subcluster.
    xt::xarray<double> a = xt::linalg::dot(this->centroids_, subcluster->centroid_);
    xt::xarray<float> dist_matrix = 1 - a / (xt::sum(this->centroids_, 1) + xt::sum(subcluster->centroid_) - a);
    int closest_index = xt::argmin(dist_matrix)();
    _CFSubcluster* closest_subcluster = this->subclusters_[closest_index];

    // If the subcluster has a child, we need a recursive strategy.
    if (closest_subcluster->child_ != nullptr) {
        bool split_child = closest_subcluster->child_->insert_cf_subcluster(subcluster);

        if (!split_child) {
            // If it is determined that the childneed not be split, we
            // can just update the closest_subcluster
            closest_subcluster->update(subcluster);
            xt::row(this->init_centroids_, closest_index) = this->subclusters_[closest_index]->centroid_;
            this->centroids_ = xt::view(this->init_centroids_, xt::range(0, this->n_samples + 1), xt::all());
            return false;
        }

        // things not too good. we need to redistribute the subclusters in
        // our child node, and add a new subcluster in the parent
        // subcluster to accommodate the new child.
        else {
            std::pair<_CFSubcluster*, _CFSubcluster*> subcluster_result = _split_node(
                closest_subcluster->child_, 
                threshold, 
                branching_factor
                );
            _CFSubcluster* new_subcluster1 = subcluster_result.first;
            _CFSubcluster* new_subcluster2 = subcluster_result.second;
            this->update_split_subclusters(
                closest_subcluster,
                new_subcluster1,
                new_subcluster2
            );

            if (this->subclusters_.size() > this->branching_factor) {
                return true;
            }
            return false;
        }
    }

    // good to go!
    else {
        bool merged = closest_subcluster->merge_subcluster(subcluster, this->threshold);
        if (merged) {
            xt::row(this->init_centroids_, closest_index) = closest_subcluster->centroid_;
            this->centroids_ = xt::view(this->init_centroids_, xt::range(0, this->n_samples + 1), xt::all());
            return false;
        }

        // not close to any other subclusters, and we still
        // have space, so add.
        else if (this->subclusters_.size() < this->branching_factor) {
            this->append_subcluster(subcluster);
            return false;
        }

        // We do not have enough space nor is it closer to an
        // other subcluster. We need to split
        else {
            this->append_subcluster(subcluster);
            return true;
        }
    }
}