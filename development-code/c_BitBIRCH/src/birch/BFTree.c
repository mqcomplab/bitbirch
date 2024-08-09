#include <stdbool.h>
#include <stdio.h>
#include "../../include/birch/BFTree.h"
#include "../../include/birch/BFNode.h"
#include "../../include/util/Smem.h"
#include "../../include/util/Integer.h"


void tree_split_root(BFTree *tree);


BFTree* tree_create
(
    int dimensionality,
    int branching_factor,
    double threshold,
    double (*distance)(struct subcluster*, struct subcluster*),
    bool apply_merging_refinement
)
{
    BFTree* tree = (BFTree*) smalloc(sizeof(BFTree));

    tree->dimensionality = dimensionality;
    tree->root = node_create(branching_factor, threshold, distance, true, apply_merging_refinement);
    tree->leaf_list = node_create(branching_factor, threshold, distance, true, apply_merging_refinement);
    tree->leaf_list->next_leaf = tree->root;
    tree->instance_index = 0;

    return tree;
}


void tree_free_rec(BFNode* root)
{
    int i;
    int subclusters_size;
    BFSubcluster *curr_subcluster;

    subclusters_size = array_size(root->subclusters);

    for (i = 0; i < subclusters_size; ++i)
    {
        curr_subcluster = (BFSubcluster*) array_get(root->subclusters, i);
        if (curr_subcluster->child != NULL)
        {
            tree_free_rec(curr_subcluster->child);
        }

        if (curr_subcluster->indexes != NULL)
        {
            array_deep_clear(curr_subcluster->indexes);
        }
        subcluster_free(curr_subcluster);
    }

    node_free(root);
}


void tree_free(BFTree* tree)
{
    tree_free_rec(tree->root);
    node_free(tree->leaf_list);
    free(tree);
}


int tree_insert(BFTree* tree, double* sample)
{
    int instance_index = tree->instance_index;
    tree->instance_index++;

    BFSubcluster* subcluster = subcluster_create(sample, tree->dimensionality, instance_index);

    bool hold_memory = false;
    bool dont_split = node_insert_subcluster(tree->root, subcluster, &hold_memory);
   
    if (dont_split == false)
    {
        // if dontSplit is false, it means there was not enough space to insert the new subcluster in the tree,
        // therefore wee need to split the root to make more room
        tree_split_root(tree);
    }


    if (hold_memory == false)
    {
        subcluster_free(subcluster);
    }
    
    return instance_index;
}


void tree_split_root(BFTree *tree)
{
    // the split happens by finding the two subclusters in this node that are the most far apart
    // we then use these two subclusters as a "pivot" to redistribute the old subclusters into two new nodes

    Pair* pair;
    BFSubcluster* new_subcluster_1;
    BFSubcluster* new_subcluster_2;
    BFNode* new_node_1;
    BFNode* new_node_2;
    BFNode* new_root;

    pair = pair_find_farthest(tree->root->subclusters, tree->root->distance);

    new_subcluster_1 = subcluster_create_default(pair->s1->dim);
    new_node_1 = node_create(tree->root->branching_factor,
                             tree->root->threshold,
                             tree->root->distance,
                             tree->root->is_leaf,
                             tree->root->apply_merging_refinement);
    new_subcluster_1->child = new_node_1;

    new_subcluster_2 = subcluster_create_default(pair->s1->dim);
    new_node_2 = node_create(tree->root->branching_factor,
                             tree->root->threshold, tree->root->distance,
                             tree->root->is_leaf,
                             tree->root->apply_merging_refinement);
    new_subcluster_2->child = new_node_2;

    // the new root that hosts the new subclusters
    new_root = node_create(tree->root->branching_factor,
                           tree->root->threshold,
                           tree->root->distance,
                           false,
                           tree->root->apply_merging_refinement);
    array_add(new_root->subclusters, new_subcluster_1);
    array_add(new_root->subclusters, new_subcluster_2);

    // this updates the pointers to the list of leaves
    if(tree->root->is_leaf == true)
    {
        // if root was a leaf
        tree->leaf_list->next_leaf = new_node_1;
        new_node_1->prev_leaf = tree->leaf_list;
        new_node_1->next_leaf = new_node_2;
        new_node_2->prev_leaf = new_node_1;
    }

    // redistributes the subclusters in the root between newsubcluster1 and newsubcluster2
    // according to the distance to p.e1 and p.e2
    node_redistribute_subclusters(tree->root, tree->root->subclusters, pair, new_subcluster_1, new_subcluster_2);

    // updates the root
    free(pair);
    node_free(tree->root);
    tree->root = new_root;
}


Array* tree_get_subclusters(BFTree* tree)
{
    Array* subclusters = array_create(1);
    BFNode* leaf = tree->leaf_list->next_leaf; // the first leaf is dummy!

    while(leaf != NULL)
    {
        if(!node_is_dummy(leaf))
        {

            for (int i = 0; i < array_size(leaf->subclusters); ++i)
            {
                BFSubcluster* subcluster = (BFSubcluster*) array_get(leaf->subclusters, i);
                array_add(subclusters, array_clone(subcluster->indexes));
            }
        }
        leaf = leaf->next_leaf;
    }

    return subclusters;
}


int tree_count_subclusters(BFTree* tree)
{
    int count = 0;
    BFNode* leaf = tree->leaf_list->next_leaf; // the first leaf is dummy!

    while(leaf != NULL)
    {
        if(!node_is_dummy(leaf))
        {
            for (int i = 0; i < array_size(leaf->subclusters); ++i)
            {
                BFSubcluster* subcluster = (BFSubcluster*) array_get(leaf->subclusters, i);
                for (int j = 0; j < array_size(subcluster->indexes); ++j)
                {
                    ++count;
                }
            }
        }
        leaf = leaf->next_leaf;
    }

    return count;
}


int* tree_get_cluster_id_by_instance_index(BFTree* tree)
{
    BFNode* leaf = tree->leaf_list->next_leaf; // the first leaf is dummy!
    int* cluster_id_by_subcluster_index = smalloc(tree->instance_index * sizeof(int));
    int cluster_id = 0;

    while(leaf != NULL)
    {
        if(!node_is_dummy(leaf))
        {
            for (int i = 0; i < array_size(leaf->subclusters); ++i)
            {
                BFSubcluster* subcluster = (BFSubcluster*) array_get(leaf->subclusters, i);
                for (int j = 0; j < array_size(subcluster->indexes); ++j)
                {
                    Integer* index = (Integer*) array_get(subcluster->indexes, j);
                    cluster_id_by_subcluster_index[index->value] = cluster_id;
                }
                ++cluster_id;
            }
        }
        leaf = leaf->next_leaf;
    }

    return cluster_id_by_subcluster_index;
}
