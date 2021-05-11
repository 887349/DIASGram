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

void matrix_init(TYPE **&data, float v = 0.0, int d = 0, int r = 0, int c = 0)
{

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                data[k][i * c + j] = rand()%10;
            }
        }
    }
}

/* copy constructor test
TYPE **copy_constr(TYPE **&that, int d = 0, int r = 0, int c = 0)
{
    TYPE **data;

    data = static_matrix_create(d, r, c);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                data[k][i * 2 + j] = that[k][i * 2 + j];
            }
        }
    }

    return data;
}
*/

/*main to test various functions*/
int main() {
    float **that, **data;

    that = static_matrix_create(3, 2, 2);
    matrix_init(that, 10, 3, 2, 2);
    data = copy_constr(that, 3, 2, 2);

    matrix_print(that, 3, 2, 2);
    std::cout << "----------" << std::endl;
    matrix_print(data, 4, 2, 2);

    matrix_del(that, 3);
    matrix_del(data, 4);

    return 0;
}