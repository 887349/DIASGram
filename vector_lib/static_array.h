#include <stdio.h>
#include <stdlib.h>
#define TYPE float

TYPE *static_create(int dim);

TYPE **static_pointer_create(int dim);

TYPE *static_matrix_create(int d, int r, int c, float v);

void v_random_filler(TYPE * v, int dim);

void static_print_union(TYPE arr1[], TYPE arr2[], int len1, int len2);

int static_print_intersection(TYPE * arr1, TYPE * arr2, int len1, int len2);

void del_at_pos(TYPE * v, int dim, int pos);

void modify_at_pos(TYPE * v, int dim, int pos, TYPE value);

void del_value(TYPE * v, int dim, TYPE value);

void modify_value(TYPE * v, int dim, TYPE previous, TYPE value);

TYPE *invert_v(TYPE * v, int dim);

void static_del(TYPE * v, int d);

int get_value_pos(TYPE * v, int dim, TYPE value);

void v_print(TYPE * v, int dim);