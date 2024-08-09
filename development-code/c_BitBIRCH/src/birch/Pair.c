#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "../../include/birch/Pair.h"
#include "../../include/util/Smem.h"


Pair* pair_create_default()
{
    return (Pair*) smalloc(sizeof(Pair));
}


Pair* pair_create(BFSubcluster* s1, BFSubcluster* s2)
{
    Pair* pair = pair_create_default();
    pair->s1 = s1;
    pair->s2 = s2;
    return pair;
}


bool pair_cmp(Pair* p1, Pair* p2)
{
    if (subcluster_cmp(p1->s1, p2->s1) == true
        && subcluster_cmp(p1->s2, p2->s2) == true)
    {
        return true;
    }

    if (subcluster_cmp(p1->s1, p2->s2) == true
        && subcluster_cmp(p1->s2, p2->s1) == true)
    {
        return true;
    }

    return false;
}

//max_separation 
/*
Pair* pair_find_farthest
(
    //each entry is a subcluster
    Array* entries,
    double (*distance)(struct entry*, struct entry*)
)
{
    int i, j;
    double max_dist;
    double curr_dist;
    Entry* e1;
    Entry* e2;
    Pair* pair;

    if (array_size(entries) < 2)
    {
        return NULL;
    }

    pair = pair_create_default();
    max_dist = -1;

    for (i = 0; i < array_size(entries) - 1; ++i)
    {
        for (j = i + 1; j < array_size(entries); ++j)
        {
            e1 = (Entry*) array_get(entries, i);
            e2 = (Entry*) array_get(entries, j);

            curr_dist = distance(e1, e2);

            if(curr_dist > max_dist)
            {
                pair->e1 = e1;
                pair->e2 = e2;
                max_dist = curr_dist;
            }
        }
    }

    return pair;
}
*/
 
//n_samp is height; n_feat is width
int *dot(int **centroids, int mlc[], int n_samp, int n_feat){
    int *arr=(int*)malloc(n_samp * sizeof(int));

    //entry numbers
    for(int i=0; i<n_samp; i++){
        //entry dimensions
        arr[i]=0;
        for(int j=0; j<n_feat; j++){

            int temp=mlc[j] * centroids[i][j];
            arr[i]+=temp;

        }
    }

    return arr;
}
int argmin(float mlc[],int size){
    int min=0;
    for(int i=1; i<size; i++){
        if(mlc[i] < mlc[min]){
            min=i;
        }
    }

    return min;
}

Pair* pair_find_farthest
(
    Array* subclusters,
    double (*distance)(struct subcluster*, struct subcluster*)
)
{
    if (array_size(subclusters) < 2)
    {
        printf("here");
        return NULL;
    }
   

    BFSubcluster* t;
    t= (BFSubcluster*) array_get(subclusters,0);
    int dim_=t->dim;
    float n_samples=(float)array_size(subclusters);

    Pair* pair=pair_create_default();

    //get centroids, linear sum and pop counts
    //int centroids[n_samples][dim_];
    int** centroids = (int**)malloc(n_samples*sizeof(int*));
    if(centroids==NULL){
	printf("not allocated");
    }

    //float linear_sum[dim_];
    int* linear_sum=(int*)malloc(dim_*sizeof(int));
    //float pop_counts[n_samples];
    float* pop_counts=(float*)malloc(n_samples*sizeof(float));

    for(int i=0; i<n_samples; i++){
        centroids[i]=(int*)malloc(dim_*sizeof(int));

        BFSubcluster* temp=(BFSubcluster*) array_get(subclusters,i);
        int pop=0;

        for(int j=0; j<dim_; j++){
            centroids[i][j] = floor(temp->ls[j]/n_samples + 0.5);

            pop+=centroids[i][j];
        }

        pop_counts[i]=pop;
    }

    for(int i=0; i<dim_; i++){
        linear_sum[i]=0;
        for(int j=0; j<n_samples; j++){
            linear_sum[i]+=centroids[j][i];
        }
    }

    //calculate med centroid
    //int med[dim_];
    int* med=(int*)malloc(dim_*sizeof(int));

    for(int i=0; i<dim_; i++){
        med[i]=floor(linear_sum[i]/n_samples + 0.5);
    }

    //calculte a_med and sims_med and mol1 index
    int* a_med=dot(centroids,med,n_samples,dim_);

    int med_sum=0;
    for(int i=0; i<dim_; i++){
        med_sum+=med[i];
    }

    //int sims_med[n_samples];
    //float sims_med[n_samples];
    float* sims_med=(float*)malloc(n_samples*sizeof(float));

    for(int i=0; i<n_samples; i++){
        sims_med[i]=a_med[i] / (med_sum + pop_counts[i] - a_med[i]);
    }

    int mol1=argmin(sims_med,n_samples);

    //calculate a_mol1 and sims_mol1 and mol2
    int* a_mol1=dot(centroids,centroids[mol1],n_samples,dim_);

    //float sims_mol1[n_samples];
    float* sims_mol1=(float*)malloc(n_samples*sizeof(float));
    
    for(int i=0; i<n_samples; i++){
        sims_mol1[i]=a_mol1[i] / (pop_counts[mol1] + pop_counts[i] - a_mol1[i]);
    }

    int mol2=argmin(sims_mol1,n_samples);

    //calculate a_mol2 and sims_mol2
    int *a_mol2=dot(centroids,centroids[mol2],n_samples,dim_);

    //int sims_mol2[n_samples];
    float* sims_mol2=(float*)malloc(n_samples*sizeof(float));

    for(int i=0; i<n_samples; i++){
        sims_mol2[i]=a_mol2[i] / (pop_counts[mol2] + pop_counts[i] - a_mol2[i]);
    }

    pair->s1=(BFSubcluster*) array_get(subclusters,mol1);
    pair->s2=(BFSubcluster*) array_get(subclusters,mol2);

    free(a_med);
    free(a_mol1);
    free(a_mol2);
    free(sims_mol2);
    free(sims_mol1);
    free(med);
    free(pop_counts);
    free(linear_sum);
    for(int i=0; i<n_samples; i++){
        free(centroids[i]);
    }
    free(centroids);

    return pair;
}

Pair* pair_find_closest
(
    Array* subclusters,
    double (*distance)(struct subcluster*, struct subcluster*)
)
{
    int i, j;
    double min_dist;
    double curr_dist;
    BFSubcluster* s1;
    BFSubcluster* s2;
    Pair* pair;

    if (array_size(subclusters) < 2)
    {
        return NULL;
    }

    pair = pair_create_default();
    min_dist = DBL_MAX;

    for (i = 0; i < array_size(subclusters) - 1; ++i)
    {
        for (j = i + 1; j < array_size(subclusters); ++j)
        {
            s1 = (BFSubcluster*) array_get(subclusters, i);
            s2 = (BFSubcluster*) array_get(subclusters, j);

            curr_dist = distance(s1, s2);

            if(curr_dist < min_dist)
            {
                pair->s1 = s1;
                pair->s2 = s2;
                min_dist = curr_dist;
            }
        }
    }

    return pair;
}
