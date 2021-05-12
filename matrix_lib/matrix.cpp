#include <iostream>
#include<cstdlib>
#include "matrix.hpp"

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

void matrix_init(TYPE **&data, int d = 0, int r = 0, int c = 0, float v = 1.0)
{

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {   
                if(v == 1.0)
                    data[k][i * c + j] = rand()%10;
                else
                    data[k][i * c + j] = v;
            }
        }
    }
}

/* copy constructor test */
TYPE **copy_constr(TYPE **&data1, float v, int d = 0, int r = 0, int c = 0)
{
    TYPE **data;

    data = static_matrix_create(d, r, c);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                data[k][i * c + j] = data1[k][i * c + j] * v;
            }
        }
    }

    return data;
}
/**/

/*main to test various functions*/
int main() {
    float **that, **data;

    that = static_matrix_create(3, 2, 2);
    matrix_init(that, 3, 2, 2);
    data = copy_constr(that, -2.7, 3, 2, 2);

    matrix_print(that, 3, 2, 2);
    matrix_print(data, 3, 2, 2);

    matrix_del(that, 3);
    matrix_del(data, 3);
    return 0;
}