#ifndef PAIR_H
#define PAIR_H

#include <stdlib.h>
#include "./BFSubcluster.h"

struct pair
{
    struct subcluster* s1;
    struct subcluster* s2;
};
typedef struct pair Pair;

Pair* pair_create_default();
Pair* pair_create(BFSubcluster* s1, BFSubcluster* s2);
bool pair_cmp(Pair* p1, Pair* p2);
int *dot(int **centroids, int mlc[], int n_samp, int n_feat);
int argmin(float mlc[],int size);
Pair* pair_find_farthest(Array* entries, double (*distance)(struct subcluster*, struct subcluster*));
Pair* pair_find_closest(Array* entries, double (*distance)(struct subcluster*, struct subcluster*));

#endif
