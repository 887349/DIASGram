#include <iostream>
#include<cstdlib>
#include "matrix.h"

TYPE **static_matrix_create(int d = 0, int r = 0, int c = 0){
    TYPE **data;
    data = new TYPE*[d]{nullptr};

    for (int i = 0; i < d; i++)
        data[i] = new TYPE[r*c];

    return data;
}

void matrix_del(TYPE **&v, int d){
    for(int i = 0; i < d; i++)
		delete[] v[i];
	delete[] v;
}

void matrix_print(TYPE **&data, int d = 0, int r = 0, int c = 0){

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                std::cout << data[k][i * c + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void matrix_init(TYPE **&data, float v = 0.0){

    for (int k = 0; k < 3; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                data[k][i * 2 + j] = rand()%1000;
            }
        }
    }
}



int main() {
    float **data;

    data = static_matrix_create(3, 2, 2);

    matrix_init(data, 10);
    matrix_print(data, 3, 2, 2);

    matrix_del(data, 3);


    return 0;
}