#include <iostream>
#include<cstdlib>
#include "matrix.hpp"

TYPE **static_matrix_create(int d = 0, int r = 0, int c = 0){
    TYPE **data;
    if (d!=0) {
        data = new TYPE*[d]{nullptr};

        for (int i = 0; i < d; i++)
            data[i] = new TYPE[r*c];
    } else {
        data = new TYPE*[0];
    }
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

void matrix_init(TYPE **&data, int d = 0, int r = 0, int c = 0, float v = 0.0)
{

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {   
                data[k][i * c + j] = v;
            }
        }
    }
}

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

void clamp(TYPE **&data, int d, int r, int c, float low, float high) {
    float *temp;

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp = &data[k][i * c + j];
                if(*temp < low)
                    *temp = low;
                else if(*temp > high)
                    *temp = high;
                else    
                    continue;
            }
        }
    }
}

/*main to test various functions*/
int main() {
    float **that, **data;

    that = static_matrix_create(3, 2, 2);
    matrix_init(that, 3, 2, 2);

    matrix_print(that, 3, 2, 2);
    clamp(that, 3, 2, 2, 1.5, 5.5);    matrix_print(that, 3, 2, 2);

    matrix_del(that, 3);
    return 0;
}