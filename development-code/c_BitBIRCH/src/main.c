/* 
 Part of the tree-management code was derived from https://github.com/douglas444/birch
 Author: Douglas Monteiro <douglas.m.cavalcanti@gmail.com>

 
 BitBIRCH is an open-source clustering module based on iSIM

 BitBIRCH is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, version 3.

 BitBIRCH is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 BitBIRCH authors (C): Vicky (Vic) Jung <jungvicky@ufl.edu>
                       Ramon Alain Miranda Quintana <ramirandaq@gmail.com>, <quintana@chem.ufl.edu>

 BitBIRCH License: LGPL-3.0 https://www.gnu.org/licenses/lgpl-3.0.en.html#license-text

 Please, cite the BitBIRCH paper:
 */


#include <time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "../include/util/Smem.h"
#include "../include/util/Instance.h"
#include "../include/util/Integer.h"
#include "../include/birch/BFTree.h"

/*
double distance(Entry* e1, Entry* e2)
{
    double dist = 0;
    for (int i = 0; i < e1->dim; ++i)
    {
        double diff = (e1->ls[i] / e1->n) - (e2->ls[i] / e2->n);
        dist += diff * diff;
    }
    return sqrt(dist);
}
*/

double distance(BFSubcluster* s1, BFSubcluster* s2){
    double sim=0;

    int sum1=0;
    int sum2=0;
    double dot=0;

    for(int i=0; i<s1->dim; i++){
        int t1=floor(s1->ls[i]/s1->n +0.5);
        int t2=floor(s1->ls[i]/s1->n +0.5);

        sum1+=t1;
        sum2+=t2;
        dot+=t1*t2;
    }

    sim=dot/(sum1+sum2-dot);

    return 1-sim;
}


void print_output(BFTree* tree, Array* indexes, char *output_file_path)
{
    FILE *f = fopen(output_file_path, "w");
    if (f == NULL)
    {
        printf("ERROR: main.c/print_output(): \"Error opening file\"\n");
        exit(1);
    }

    int* cluster_id_by_subcluster_index = tree_get_cluster_id_by_instance_index(tree);

    for (int i = 0; i < array_size(indexes); ++i)
    {
        Integer* subcluster_index = (Integer*) array_get(indexes, i);
        fprintf(f, "cluster_%d\n", cluster_id_by_subcluster_index[subcluster_index->value]);
    }

    fclose(f);

    free(cluster_id_by_subcluster_index);
}


void print_arguments_message()
{
    printf("\n>>> Invalid numbers of arguments! <<<\n");
    printf("\nThe following parameters must be informed in this order and separated by whitespace:\n");
    printf("\n1. Branching Factor: An integer value;\n");
    printf("\n2. Threshold: A float value;\n");
    printf("\n3. Apply Merging Refinement: 1 for yes, 0 for no;\n");
    printf("\n4. Dataset File Path: String;\n");
    printf("\n5. Dataset Delimiters: String;\n");
    printf("\n6. Ignore Last Column of the Dataset: 1 for yes, 0 for no;\n");
    printf("\n7. Output File Name: String.\n");
    printf("\n\nThe line bellow is an example of a valid command line for running this program:\n");
    printf("./main 100 0.8 1 IRIS.csv , 1 output.csv\n\n");
}


int main(int argc, char* argv[])
{
    if (argc < 8)
    {
        print_arguments_message();
        exit(1);
    }

    int branching_factor = atoi(argv[1]);
    double threshold = atof(argv[2]);
    int apply_merging_refinement = atoi(argv[3]);
    char* input_file_path = argv[4];
    char* column_delimiter = argv[5];
    int last_column_is_label = atoi(argv[6]);
    char* output_file_path = argv[7];

    FILE* stream = fopen(input_file_path, "r");
    if (stream == NULL)
    {
        printf("ERROR: main.c/main(): \"Error opening file\"\n");
        exit(1);
    }

    char* delimiters = smalloc(sizeof(char*) * (strlen(column_delimiter) + 3));
    strcpy(delimiters, column_delimiter);
    strcat(delimiters, "\r\n");

    char line[1024];
    fgets(line, 1024, stream);
    int dimensionality = instance_calculate_dimensionality(line, delimiters, last_column_is_label);
    BFTree* tree = tree_create(dimensionality, branching_factor, threshold, distance, apply_merging_refinement == 1);
    Array* instances_indexes = array_create(1);

    do {

        double* instance = instance_read(line, dimensionality, delimiters);
        int instance_index = tree_insert(tree, instance);

        array_add(instances_indexes, integer_create(instance_index));

        free(instance);

    } while(fgets(line, 1024, stream));
    
    char buf[100];

    printf("done");
    
    fclose(stream);

    print_output(tree, instances_indexes, output_file_path);

    array_deep_clear(instances_indexes);
    array_free(instances_indexes);
    free(delimiters);
    tree_free(tree);

    return 0;
}
