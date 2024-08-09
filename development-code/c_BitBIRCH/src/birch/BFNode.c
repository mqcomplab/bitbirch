#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <stdbool.h>
#include "../../include/birch/BFSubcluster.h"
#include "../../include/birch/BFNode.h"
#include "../../include/util/Smem.h"


BFSubcluster* find_closest_subcluster(BFNode *node, BFSubcluster* subcluster);
void merging_refinement(BFNode* node, Pair* split_subclusters);
int find_closest_subcluster_(BFNode* node, BFSubcluster* subcluster);
Pair* split_subcluster(BFNode* node, BFSubcluster* closest_subcluster);
void replace_closest_pair_with_new_subclusters(BFNode *node, Pair* pair, BFSubcluster* new_subcluster_1, BFSubcluster* new_subcluster_2);
void redistribute_subclusters(BFNode* node, Array* old_subclusters_1, Array* old_subclusters_2, Pair* close_subclusters, BFSubcluster* new_subcluster_1, BFSubcluster* new_subcluster_2);
void redistribute_subclusters_merge(BFNode *node, Array* old_subclusters_1, Array* old_subclusters_2, BFSubcluster* new_subcluster);
void replace_closest_pair_with_new_merged_subcluster(BFNode* node, Pair* pair, BFSubcluster* new_subcluster);


BFNode* node_create
(
    int branching_factor,
    double threshold,
    double (*distance)(struct subcluster*, struct subcluster*),
    bool is_leaf,
    bool apply_merging_refinement
)
{
    BFNode* node = (BFNode*) smalloc(sizeof(BFNode));
    node->branching_factor = branching_factor;
    node->threshold = threshold;
    node->distance = distance;
    node->subclusters = array_create(branching_factor);
    node->is_leaf = is_leaf;
    node->next_leaf = NULL;
    node->prev_leaf = NULL;
    node->apply_merging_refinement = apply_merging_refinement;
    return node;
}


void node_free(BFNode *node)
{
    array_free(node->subclusters);
    free(node);
}


bool node_is_dummy(BFNode* node)
{
    if ((node->prev_leaf != NULL || node->next_leaf != NULL) &&
            node->branching_factor == 0 && node->threshold == 0 &&
            array_size(node->subclusters))
    {
        return true;
    }
    else
    {
        return false;
    }
}


BFSubcluster* find_closest_subcluster(BFNode *node, BFSubcluster* subcluster)
{
    int i;
    double min_dist;
    double curr_dist;
    BFSubcluster* curr_subcluster;
    BFSubcluster* closest_subcluster;

    closest_subcluster = NULL;
    min_dist = DBL_MAX;

    for (i = 0; i < array_size(node->subclusters); ++i)
    {
        curr_subcluster = (BFSubcluster*) array_get(node->subclusters, i);
        curr_dist = node->distance(curr_subcluster, subcluster);

        if (curr_dist < min_dist)
        {
            min_dist = curr_dist;
            closest_subcluster = curr_subcluster;
        }
    }

    return closest_subcluster;
}


//not used
int find_closest_subcluster_(BFNode* node, BFSubcluster* subcluster)
{
    BFSubcluster* closest_subcluster = find_closest_subcluster(node, subcluster);

    if(closest_subcluster->child == NULL)
    {
        return closest_subcluster->subcluster_id;
    }

    //return find_closest_subcluster(closest_subcluster->child, subcluster);
    return 0;
}


Pair* split_subcluster(BFNode* node, BFSubcluster* closest_subcluster)
{
    // IF there was a child, but we could not insert the new subcluster without problems THAN
    // split the child of closest subcluster

    Array* old_subclusters;
    Pair* farthest_pair;
    Pair* new_pair;
    BFNode* old_node;
    BFNode* new_node_1;
    BFNode* new_node_2;
    BFNode* prev;
    BFNode* next;
    BFSubcluster* new_subcluster_1;
    BFSubcluster* new_subcluster_2;

    old_node = closest_subcluster->child;
    old_subclusters = closest_subcluster->child->subclusters;
    farthest_pair = pair_find_farthest(old_subclusters, node->distance);

    new_subcluster_1 = subcluster_create_default(closest_subcluster->dim);
    new_node_1 = node_create(node->branching_factor, node->threshold, node->distance, old_node->is_leaf, node->apply_merging_refinement);
    new_subcluster_1->child = new_node_1;

    new_subcluster_2 = subcluster_create_default(closest_subcluster->dim);
    new_node_2 = node_create(node->branching_factor, node->threshold, node->distance, old_node->is_leaf, node->apply_merging_refinement);
    new_subcluster_2->child = new_node_2;


    if(old_node->is_leaf)
    {
        // we do this to preserve the pointers in the leafList

        prev = old_node->prev_leaf;
        next = old_node->next_leaf;

        if(prev != NULL)
        {
            prev->next_leaf = new_node_1;
        }

        if(next != NULL)
        {
            next->prev_leaf = new_node_2;
        }

        new_node_1->prev_leaf = prev;
        new_node_1->next_leaf = new_node_2;
        new_node_2->prev_leaf = new_node_1;
        new_node_2->next_leaf = next;
    }

    node_redistribute_subclusters(node, old_subclusters, farthest_pair, new_subcluster_1, new_subcluster_2);
    // redistributes the subclusters in n between newsubcluster1 and newsubcluster2
    // according to the distance to p.e1 and p.e2

    node_free(old_node);

    free(farthest_pair);

    subcluster_remove(node->subclusters, closest_subcluster); // this will be substitute by two new subclusters

    if (closest_subcluster->indexes != NULL)
    {
        array_deep_clear(closest_subcluster->indexes);
    }

    subcluster_free(closest_subcluster);

    array_add(node->subclusters, new_subcluster_1);
    array_add(node->subclusters, new_subcluster_2);

    new_pair = pair_create(new_subcluster_1, new_subcluster_2);

    return new_pair;
}


void node_redistribute_subclusters
(
    BFNode* node,
    Array* old_subclusters,
    Pair *far_subclusters,
    BFSubcluster* new_subcluster_1,
    BFSubcluster* new_subcluster_2
)
{
    int i;
    double dist_1;
    double dist_2;
    BFSubcluster* curr_subcluster;

    for (i = 0; i < array_size(old_subclusters); ++i)
    {
        curr_subcluster = (BFSubcluster*) array_get(old_subclusters, i);
        dist_1 = node->distance(far_subclusters->s1, curr_subcluster);
        dist_2 = node->distance(far_subclusters->s2, curr_subcluster);

        if (dist_1 <= dist_2)
        {
            array_add(new_subcluster_1->child->subclusters, curr_subcluster);
            subcluster_update(new_subcluster_1, curr_subcluster);
        }
        else
        {
            array_add(new_subcluster_2->child->subclusters, curr_subcluster);
            subcluster_update(new_subcluster_2, curr_subcluster);
        }
    }
}


//inferior one
void redistribute_subclusters
(
    BFNode* node,
    Array* old_subclusters_1,
    Array* old_subclusters_2,
    Pair* close_subclusters,
    BFSubcluster* new_subcluster_1,
    BFSubcluster* new_subcluster_2
)
{
    int i;
    Array* v;
    BFSubcluster* curr_subcluster;
    double dist_1;
    double dist_2;

    v = array_create(node->branching_factor * 2);
    array_add_all(v, old_subclusters_1);
    array_add_all(v, old_subclusters_2);

    for (i = 0; i < array_size(v); ++i)
    {
        curr_subcluster = (BFSubcluster*) array_get(v, i);
        dist_1 = node->distance(close_subclusters->s1, curr_subcluster);
        dist_2 = node->distance(close_subclusters->s2, curr_subcluster);

        if (dist_1 <= dist_2)
        {
            if(array_size(new_subcluster_1->child->subclusters) < node->branching_factor)
            {
                array_add(new_subcluster_1->child->subclusters, curr_subcluster);
                subcluster_update(new_subcluster_1, curr_subcluster);
            }
            else
            {
                array_add(new_subcluster_2->child->subclusters, curr_subcluster);
                subcluster_update(new_subcluster_2, curr_subcluster);
            }
        }
        else if(dist_2 < dist_1)
        {
            if(array_size(new_subcluster_2->child->subclusters) < node->branching_factor)
            {
                array_add(new_subcluster_2->child->subclusters, curr_subcluster);
                subcluster_update(new_subcluster_2, curr_subcluster);
            }
            else
            {
                array_add(new_subcluster_1->child->subclusters, curr_subcluster);
                subcluster_update(new_subcluster_1, curr_subcluster);
            }
        }
    }
}


//only fore merging refinement
void redistribute_subclusters_merge
(
    BFNode *node,
    Array* old_subclusters_1,
    Array* old_subclusters_2,
    BFSubcluster* new_subcluster
)
{
    int i;
    Array* v;
    BFSubcluster* curr_subcluster;

    v = array_create(node->branching_factor);
    array_add_all(v, old_subclusters_1);
    array_add_all(v, old_subclusters_2);

    for (i = 0; i < array_size(v); ++i)
    {
        curr_subcluster = (BFSubcluster*) array_get(v, i);
        array_add(new_subcluster->child->subclusters, curr_subcluster);
        subcluster_update(new_subcluster, curr_subcluster);
    }
}


//only for merging refinement
void replace_closest_pair_with_new_subclusters
(
    BFNode *node,
    Pair* pair,
    BFSubcluster* new_subcluster_1,
    BFSubcluster* new_subcluster_2
)
{
    int i;

    for (i = 0; i < array_size(node->subclusters); i++)
    {
        if(subcluster_cmp((BFSubcluster*) array_get(node->subclusters, i), pair->s1) == true)
        {
            array_set(node->subclusters, i, new_subcluster_1);
        }
        else if(subcluster_cmp((BFSubcluster*) array_get(node->subclusters, i), pair->s2) == true)
        {
            array_set(node->subclusters, i, new_subcluster_2);
        }
    }
}


//only for merging refinement
void replace_closest_pair_with_new_merged_subcluster
(
    BFNode* node,
    Pair* pair,
    BFSubcluster* new_subcluster
)
{
    int i;

    for (i = 0; i < array_size(node->subclusters); i++)
    {
        if (subcluster_cmp((BFSubcluster*) array_get(node->subclusters, i), pair->s1) == true)
        {
            array_set(node->subclusters, i, new_subcluster);
        }
        else if (subcluster_cmp((BFSubcluster*) array_get(node->subclusters, i), pair->s2) == true)
        {
            array_remove_by_index(node->subclusters, i);
        }
    }
}


void merging_refinement(BFNode* node, Pair* split_subclusters)
{

    BFNode* old_node_1;
    BFNode* old_node_2;
    BFNode* new_node_1;
    BFNode* new_node_2;
    BFNode* new_node;
    BFNode* dummy_node;
    Array* old_node_1_subclusters;
    Array* old_node_2_subclusters;
    Array* node_subclusters;
    BFSubcluster* new_subcluster_1;
    BFSubcluster* new_subcluster_2;
    BFSubcluster* new_subcluster;
    Pair* pair;

    node_subclusters = node->subclusters;
    pair = pair_find_closest(node_subclusters, node->distance);

    if (pair == NULL)
    {
        // not possible to find a valid pair
        return;
    }

    if (pair_cmp(pair, split_subclusters))
    {
        return; // if the closet pair is the one that was just split, we terminate
    }

    old_node_1 = pair->s1->child;
    old_node_2 = pair->s2->child;

    old_node_1_subclusters = array_clone(old_node_1->subclusters);
    old_node_2_subclusters = array_clone(old_node_2->subclusters);

    if(old_node_1->is_leaf != old_node_2->is_leaf)
    {
        // just to make sure everything is going ok
        printf("ERROR: node.c/merging_refinement(): \"Nodes at the same level must have same leaf status\"\n");
        exit(1);
    }

    if((array_size(old_node_1_subclusters) + array_size(old_node_2_subclusters)) > node->branching_factor)
    {
        // the two nodes cannot be merged into one (they will not fit)
        // in this case we simply redistribute them between p.e1 and p.e2

        new_subcluster_1 = subcluster_create_default(split_subclusters->s1->dim);
        // note: in the CFNode construction below the last parameter is false
        // because a split cannot happen at the leaf level
        // (the only exception is when the root is first split, but that's treated separately)
        new_node_1 = old_node_1;
        array_clear(new_node_1->subclusters);
        new_subcluster_1->child = new_node_1;

        new_subcluster_2 = subcluster_create_default(split_subclusters->s1->dim);
        new_node_2 = old_node_2;
        array_clear(new_node_2->subclusters);
        new_subcluster_2->child = new_node_2;

        redistribute_subclusters(node, old_node_1_subclusters, old_node_2_subclusters, pair, new_subcluster_1, new_subcluster_2);
        replace_closest_pair_with_new_subclusters(node, pair, new_subcluster_1, new_subcluster_2);
    }
    else
    {
        // if the the two closest subclusters can actually be merged into one single subcluster

        new_subcluster = subcluster_create_default(split_subclusters->s1->dim);
        // note: in the CFNode construction below the last parameter is false
        // because a split cannot happen at the leaf level
        // (the only exception is when the root is first split, but that's treated separately)
        new_node = node_create(node->branching_factor, node->threshold, node->distance, old_node_1->is_leaf, node->apply_merging_refinement);
        new_subcluster->child = new_node;

        redistribute_subclusters_merge(node, old_node_1_subclusters, old_node_2_subclusters, new_subcluster);

        if (old_node_1->is_leaf && old_node_2->is_leaf)
        {
            // this is done to maintain proper links in the leafList
            if(old_node_1->prev_leaf != NULL)
            {
                old_node_1->prev_leaf->next_leaf = new_node;
            }
            if(old_node_1->next_leaf != NULL)
            {
                old_node_1->next_leaf->prev_leaf = new_node;
            }
            new_node->prev_leaf = old_node_1->prev_leaf;
            new_node->next_leaf = old_node_1->next_leaf;

            // this is a dummy node that is only used to maintain proper links in the leafList
            // no CFsubcluster will ever point to this leaf
            dummy_node = node_create(0, 0, NULL, true, false);

            if (old_node_2->prev_leaf != NULL)
            {
                old_node_2->prev_leaf->next_leaf = dummy_node;
            }
            if (old_node_2->next_leaf != NULL)
            {
                old_node_2->next_leaf->prev_leaf = dummy_node;
            }

            dummy_node->prev_leaf = old_node_2->prev_leaf;
            dummy_node->next_leaf = old_node_2->next_leaf;
        }

        replace_closest_pair_with_new_merged_subcluster(node, pair, new_subcluster);
    }

    // merging refinement is done
}


bool node_insert_subcluster(BFNode* node, BFSubcluster* subcluster, bool* hold_memory)
{

    BFSubcluster* closest_subcluster;
    bool dont_split;
    Pair* split_pair;

    if(array_size(node->subclusters) == 0)
    {
        *hold_memory = true;
        array_add(node->subclusters, subcluster);
        return true;
    }

    closest_subcluster = find_closest_subcluster(node, subcluster);
    
    /*for(int i=0; i<closest_subcluster->dim;i++){
        printf("%d",closest_subcluster->ls[i]);
    }
    printf("\n");
    */

    dont_split = false;

    if(closest_subcluster->child != NULL)
    {
        dont_split = node_insert_subcluster(closest_subcluster->child, subcluster, hold_memory);

        if(dont_split == true)
        {
            subcluster_update(closest_subcluster, subcluster); // this updates the CF to reflect the additional subcluster
            return true;
        }
        else
        {
            // if the node below /closest/ didn't have enough room to host the new subcluster
            // we need to split it
            split_pair = split_subcluster(node, closest_subcluster);

            // after adding the new subclusters derived from splitting /closest/ to this node,
            // if we have more than maxsubclusters we return false,
            // so that the parent node will be split as well to redistribute the "load"
            if(array_size(node->subclusters) > node->branching_factor)
            {
                free(split_pair);
                return false;
            }
            else
            {
                // splitting stops at this node
                if(node->apply_merging_refinement)
                {
                    // performs step 4 of insert process (see BIRCH paper, Section 4.3)
                    merging_refinement(node, split_pair);
                }
                free(split_pair);
                return true;
            }
        }
    }
    else if(subcluster_is_within_threshold(closest_subcluster, subcluster, node->threshold, node->distance))
    {
        // if  dist(closest,e) <= T, /e/ will be "absorbed" by /closest/
        subcluster_update(closest_subcluster, subcluster);
        return true; // no split necessary at the parent level
    }
    else if(array_size(node->subclusters) < node->branching_factor)
    {
        // if /closest/ does not have children, and dist(closest,e) > T
        // if there is enough room in this node, we simply add e to it
        *hold_memory = true;
        array_add(node->subclusters, subcluster);
        return true; // no split necessary at the parent level
    }
    else
    {
        // not enough space on this node
        *hold_memory = true;
        array_add(node->subclusters, subcluster); // adds it momentarily to this node
        return false;   // returns false so that the parent subcluster will be split
    }

}
