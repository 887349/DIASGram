#include <stdio.h>
#include <stdlib.h>
#include <string>
#include<cstdlib>
#define TYPE float

using namespace std;

TYPE **static_matrix_create(int d, int r, int c);

void matrix_del(TYPE **&data, int d);

void matrixc_print(TYPE **&data, int d, int r, int c);

void matrix_init(TYPE **&data, int d, int r, int c, float v);

TYPE **copy_constr(TYPE **&data1, float v, int d, int r, int c);

void write_file(string filename, TYPE **&data, int d, int r, int c);

float getMax(TYPE **&data, int k, int r, int c, int d);

float getMax(TYPE **&data, int k, int r, int c, int d);

void swap_channel(TYPE **&data, int a, int b);