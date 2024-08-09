#ifndef NODE_H
#define NODE_H

#include <stdlib.h>
#include "./Pair.h"
#include "../util/Array.h"

struct node
{
    Array *subclusters;
    int branching_factor;
    double threshold;
    bool is_leaf;
    double (*distance)(struct subcluster*, struct subcluster*);
    struct node *next_leaf;
    struct node *prev_leaf;
    bool apply_merging_refinement;
};
typedef struct node BFNode;

void node_free(BFNode *node);
BFNode* node_create(int branching_factor, double threshold, double (*distance)(struct subcluster*, struct subcluster*), bool is_leaf, bool apply_merging_refinement);
bool node_is_dummy(BFNode* node);
void node_redistribute_subclusters(BFNode* node, Array* old_subclusters, Pair *far_subclusters, BFSubcluster* new_subcluster_1, BFSubcluster* new_subcluster_2);
bool node_insert_subcluster(BFNode* node, BFSubcluster* subcluster, bool* hold_memory);

#endif
