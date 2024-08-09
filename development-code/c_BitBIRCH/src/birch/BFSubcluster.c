#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/birch/BFSubcluster.h"
#include "../../include/util/Smem.h"
#include "../../include/util/Integer.h"
#include "../../include/util/Util.h"


#define INDEXES_INITIAL_SIZE 4


BFSubcluster* subcluster_create_default(int dim)
{
    BFSubcluster* subcluster = (BFSubcluster*) smalloc(sizeof(BFSubcluster));

    subcluster->dim = dim;
    subcluster->n = 0;
    subcluster->ls = (double*) scalloc(dim, sizeof(double));
    subcluster->ss = (double*) scalloc(dim, sizeof(double));
    subcluster->child = NULL;
    subcluster->indexes = NULL;
    subcluster->subcluster_id = -1;

    return subcluster;
}


BFSubcluster* subcluster_create
(
    double* x,
    int dim,
    int index
)
{
    int i;

    BFSubcluster* subcluster = subcluster_create_default(dim);

    subcluster->n = 1;
    subcluster->dim = dim;

    smemcpy(subcluster->ls, x, dim * sizeof(x));

    for (i = 0; i < subcluster->dim; i++)
    {
        subcluster->ss[i] = subcluster->ls[i] * subcluster->ls[i];
    }

    subcluster->indexes = array_create(INDEXES_INITIAL_SIZE);
    array_add(subcluster->indexes, integer_create(index));

    return subcluster;

}


void subcluster_free(BFSubcluster *subcluster)
{
    subcluster->dim = 0;
    subcluster->n = 0;

    if (subcluster->ls != NULL)
    {
        free(subcluster->ls);
    }

    if (subcluster->ss != NULL)
    {
        free(subcluster->ss);
    }

    if (subcluster->indexes != NULL)
    {
        array_free(subcluster->indexes);
    }

    free(subcluster);

}


void subcluster_update(BFSubcluster *s1, BFSubcluster *s2)
{
    int i;

    s1->n += s2->n;

    for(i = 0; i < s1->dim; ++i)
    {
        s1->ls[i] += s2->ls[i];
    }

    for(i = 0; i < s1->dim; i++)
    {
        s1->ss[i] += s2->ss[i];
    }

    if(s1->child == NULL)
    {
        if(s1->indexes != NULL && s2->indexes != NULL)
        {
            array_add_all(s1->indexes, s2->indexes);
        }
        else if(s1->indexes == NULL && s2->indexes != NULL)
        {
            s1->indexes = array_clone(s2->indexes);
        }
    }


}


//here
/*
bool entry_is_within_threshold
(
    Entry *e1,
    Entry *e2,
    double threshold,
    double (*distance)(Entry*, Entry*)
)
{
    double dist = distance(e1, e2);

    if(dist == 0 || dist <= threshold)
    {
        return true;
    }
    else
    {
        return false;
    }
}
*/

double jt_isim(int s1[], int n, int size)
{
    int sum_kq=0;
    int sum_kqsq=0; 

    for(int i=0; i<size; i++){
        sum_kq+=s1[i];
        sum_kqsq+=s1[i] * s1[i];
    }

    double a=(sum_kqsq-sum_kq)/2.0;

    return a/(a+n*sum_kq-sum_kqsq);
}
bool subcluster_is_within_threshold
(
    BFSubcluster *s1,
    BFSubcluster *s2,
    double threshold,
    double(*distance)(BFSubcluster*,BFSubcluster*)
)
{
    int size=s1->dim;
    int a1[size];
    int a2[size];
    int new_n=s1->n + s2->n;

    for(int i=0; i<size; i++){
        int new_ls = s1->ls[i] + s2->ls[i];
        int new_cent=floor(new_ls/new_n + 0.5);

        a1[i]=new_ls + new_cent;
        a2[i] = new_ls;
    }

    double j1=jt_isim(a1,new_n+1,size) * (new_n+1);
    double j2=jt_isim(a2,new_n,size) * (new_n-1);

    if(j1-j2 >= (threshold*2)){
        return true;
    }
    return false;
}

bool subcluster_cmp(BFSubcluster *s1, BFSubcluster *s2)
{
    if (s1->n != s2->n)
        return false;

    if (s1->child != NULL && s2->child == NULL)
        return false;

    if (s1->child == NULL && s2->child != NULL)
        return false;

    if(s1->child != NULL && s1->child != s2->child)
        return false;

    if (s1->indexes == NULL && s2->indexes != NULL)
        return false;

    if (s1->indexes != NULL && s2->indexes == NULL)
        return false;

    if (double_cmp(s1->ls, s2->ls, s1->dim) == false)
        return false;

    if (double_cmp(s1->ss, s2->ss, s1->dim) == false)
        return false;

    if(s1->indexes != NULL && integer_array_cmp(s1->indexes, s2->indexes) == false)
        return false;

    return true;
}


void subcluster_remove(Array* subclusters, BFSubcluster* subcluster)
{
    int index = 0;
    while (index < array_size(subclusters) && subcluster_cmp(subcluster, array_get(subclusters, index)) == false)
    {
        ++index;
    }

    if (index < array_size(subclusters))
    {
        array_remove_by_index(subclusters, index);
    }
}
