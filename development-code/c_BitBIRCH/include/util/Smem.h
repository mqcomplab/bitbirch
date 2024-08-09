#ifndef SMEM_H
#define SMEM_H

#include <stdlib.h>

void* smalloc(size_t mem_size);
void* scalloc(int num, size_t elem_size);
void* srealloc(void *mem_pos, size_t mem_size);
void* smemcpy(void *mem_pos_dest, void *mem_pos_src, size_t mem_size);

#endif

