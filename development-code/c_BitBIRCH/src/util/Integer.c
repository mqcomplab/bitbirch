#include "../../include/util/Integer.h"
#include "../../include/util/Smem.h"


Integer* integer_create(int value)
{
    Integer* integer = (Integer*) smalloc(sizeof(Integer));
    integer->value = value;
    return integer;
}


bool integer_array_cmp(Array* a1, Array* a2)
{
    for (int i = 0; i < array_size(a1); ++i)
    {
        if (((Integer*) array_get(a1, i))->value
            != ((Integer*) array_get(a2, i))->value)
        {
            return false;
        }
    }
    return true;
}
