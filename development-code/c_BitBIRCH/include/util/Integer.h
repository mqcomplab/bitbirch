#ifndef INTEGER_H
#define INTEGER_H

#include <stdbool.h>
#include "Array.h"

struct integer {
    int value;
};
typedef struct integer Integer;

Integer* integer_create(int value);
bool integer_array_cmp(Array* a1, Array* a2);

#endif
