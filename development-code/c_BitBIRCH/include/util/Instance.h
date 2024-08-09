#ifndef SAMPLE_H
#define SAMPLE_H

#include <stdbool.h>
#include "Array.h"

double* instance_read(char* line, int dimensionality, char* delimiters);
int instance_calculate_dimensionality(char* line, char* delimiters, bool last_column_is_label);

#endif

