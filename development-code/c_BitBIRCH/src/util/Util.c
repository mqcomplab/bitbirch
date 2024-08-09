#include <stdbool.h>
#include "../../include/util/Util.h"

bool double_cmp(double* d1, double* d2, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (d1[i] != d2[i])
        {
            return false;
        }
    }
    return true;
}
