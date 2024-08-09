#ifndef TREE_H
#define TREE_H

#include <stdlib.h>
#include "BFSubcluster.h"

struct tree
{
    struct node* root;
    struct node* leaf_list;
    int instance_index;
    int dimensionality;
};
typedef struct tree BFTree;

int* tree_get_cluster_id_by_instance_index(BFTree* tree);
BFTree* tree_create(int dimensionality, int branching_factor, double threshold, double (*distance)(BFSubcluster*, BFSubcluster*), bool apply_merging_refinement);
int tree_insert(BFTree* tree, double* sample);
void tree_free(BFTree* tree);

#endif
