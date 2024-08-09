#ifndef ENTRY_H
#define ENTRY_H

#include "../util/Array.h"

struct subcluster
{
    int dim;
    int n;
    double* ls;
    double* ss;
    struct node* child;
    Array* indexes;
    int subcluster_id;
};
typedef struct subcluster BFSubcluster;

BFSubcluster* subcluster_create_default(int dim);
BFSubcluster* subcluster_create(double* x, int dim, int index);
void subcluster_free(BFSubcluster *subcluster);
void subcluster_update(BFSubcluster *s1, BFSubcluster *s2);
double jt_isim(int s1[], int n, int size);
bool subcluster_is_within_threshold(BFSubcluster *s1, BFSubcluster *s2, double threshold, double (*distance)(BFSubcluster*, BFSubcluster*));
bool subcluster_cmp(BFSubcluster *s1, BFSubcluster *s2);
void subcluster_remove(Array* subclusters, BFSubcluster* subcluster);

#endif
