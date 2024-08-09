#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "../../include/util/Smem.h"
#include "../../include/util/Instance.h"

char* get_field(char* line, char* delimiters, int num)
{
    char* tok;
    char* tmp = strdup(line);
    for (tok = strtok(tmp, delimiters); tok && *tok; tok = strtok(NULL, delimiters))
    {
        if (!((--num) + 1))
        {
            char* field = strdup(tok);
            //printf("%c",*tok);
            free(tmp);
            return field;
        }
    }
    free(tmp);
    return NULL;
}

double* instance_read(char* line, int dimensionality, char* delimiters)
{
    double* sample = smalloc(sizeof(double*) * dimensionality);

    for (int i = 0; i < dimensionality; ++i)
    {
        char* field = get_field(line, delimiters, i);
        sscanf(field, "%lf", &sample[i]);
        free(field);
    }
    
    return sample;
}

int instance_calculate_dimensionality(char* line, char* delimiters, bool last_column_is_label)
{
    char *p;
    int count = 0;
    char* tmp = strdup(line);

    p = strtok(tmp, delimiters);

    while (p != NULL)
    {
        count++;
        p = strtok(NULL, delimiters);
    }

    free(tmp);

    return (last_column_is_label == true) ? count - 1 : count;
}
