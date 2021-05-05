#include <stdio.h>
#include <stdlib.h>
#define TYPE int

typedef struct dynamic_array{
    TYPE *values;
    unsigned int available;
    unsigned int used;
}dynamic_array_t;

dynamic_array_t *dynamic_create(int dim);

void dynamic_random_filler(dynamic_array_t *v);

void dynamic_print(dynamic_array_t *v);

void dynamic_del(dynamic_array_t *v);

void dynamic_sharp_amp(dynamic_array_t *v);

void dynamic_amp(dynamic_array_t *v);

void dynamic_bulk_amp(dynamic_array_t *v);

void dynamic_add(dynamic_array_t *v, TYPE value);

void dynamic_modify_value(dynamic_array_t *v, int previous, int value);

void dynamic_modify_at_position(dynamic_array_t *v, unsigned int pos, TYPE value);

unsigned int dynamic_get_position(dynamic_array_t *v, TYPE value);

void dynamic_print_intersection(dynamic_array_t *arr1, dynamic_array_t *arr2);

void dynamic_print_union(dynamic_array_t *arr1, dynamic_array_t *arr2);

void dynamic_del_value(dynamic_array_t *v, TYPE value);