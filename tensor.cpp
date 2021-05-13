#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"
#include "matrix_lib/matrix.hpp"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;


/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}
 
Tensor::Tensor(){
    data = static_matrix_create(0, 0, 0);

    r = 0; // number of rows
    c = 0; // number of columns
    d = 0; // tensor depth
}

Tensor::Tensor(int r, int c, int d, float v=0.0) {
    data = static_matrix_create(d, r, c);
    matrix_init(data, v);

    r = r; // number of rows
    c = c; // number of columns
    d = d; // tensor depth
}

Tensor::Tensor(const Tensor& that) {

    data = static_matrix_create(that.d, that.r, that.c);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                data[k][i * c + j] = that.data[k][i * c + j];
            }
        }
    }

    r = that.r;
    c = that.c;
    d = that.d;
}

Tensor::~Tensor() {
    matrix_del(data, d);
}

float Tensor::operator()(int i, int j, int k) const {
    float res{0.0};

    if(k < d || k < d)
        throw(index_out_of_bound());
    else if (i < r || i < r)
        throw(index_out_of_bound());
    else if (j < c || j < c)
        throw(index_out_of_bound());
    else
        res = data[k][(i*c)+j];

    return res;
}

float& Tensor::operator()(int i, int j, int k) {
    float *res;

    if (i < d || i < d)
        throw(index_out_of_bound());
    else if (j < r || j < r)
        throw(index_out_of_bound());
    else if (i < c || i < c)
        throw(index_out_of_bound());
    else
        return data[k][(i * c) + j];
}

Tensor Tensor::operator-(const Tensor &rhs)const {
    if (c < rhs.c || c > rhs.c)
        throw(dimension_mismatch());
    else if (r < rhs.r || r > rhs.r)
        throw(dimension_mismatch());
    else if (d < rhs.d || d > rhs.d)
        throw(dimension_mismatch());

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] = (rhs.data[k][i * c + j]) - (this->data[k][i * c + j]);
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator +(const Tensor &rhs)const {
    if (c < rhs.c || c > rhs.c)
        throw(dimension_mismatch());
    else if (r < rhs.r || r > rhs.r)
        throw(dimension_mismatch());
    else if (d < rhs.d || d > rhs.d)
        throw(dimension_mismatch());

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                temp.data[k][i * c + j] = (rhs.data[k][i * c + j]) + (this->data[k][i * c + j]);
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator*(const Tensor &rhs)const
{
    if (c < rhs.c || c > rhs.c)
        throw(dimension_mismatch());
    else if (r < rhs.r || r > rhs.r)
        throw(dimension_mismatch());
    else if (d < rhs.d || d > rhs.d)
        throw(dimension_mismatch());

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                temp.data[k][i * c + j] = (rhs.data[k][i * c + j]) * (this->data[k][i * c + j]);
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator/(const Tensor &rhs)const
{
    if (c < rhs.c || c > rhs.c)
        throw(dimension_mismatch());
    else if (r < rhs.r || r > rhs.r)
        throw(dimension_mismatch());
    else if (d < rhs.d || d > rhs.d)
        throw(dimension_mismatch());

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                temp.data[k][i * c + j] = (rhs.data[k][i * c + j]) / (this->data[k][i * c + j]);
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator-(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] - rhs;
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator+(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] + rhs;
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator*(const float &rhs)const {

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                temp.data[k][i * c + j] * rhs;
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor Tensor::operator/(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] / rhs;
            }
        }
    }

    temp.d = d;
    temp.r = r;
    temp.c = c;

    return temp;
}

Tensor &Tensor::operator=(const Tensor &other) {
    
    matrix_del(this->data, this->d);
    this->data = static_matrix_create(other.d, other.r, other.c);
    this-> d = other.d;
    this-> r = other.r;
    this-> c = other.c;

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                if (this->data[k][i * c + j] != (other.data[k][i * c + j]))
                    this->data[k][i * c + j] = (other.data[k][i * c + j]);
            }
        }
    }
    
    return *this;
}

void Tensor::init(int r, int c, int d, float v=0.0) {

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                data[k][i * c + j] * v;
            }
        }
    }
}

void Tensor::clamp(float low, float high) {
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

void Tensor::rescale(float new_max=1.0) {
    for (int k = 0; k < d; k++)
    {
        float min = getMin(k);
        float max = getMax(k);
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                if (max==min)
                    this->data[k][i * c + j] = new_max;
                else
                    this->data[k][i * c + j] = ((this->data[k][i * c + j] - min) / (max - min)) * new_max;
                
                if (this->data[k][i * c + j] < 0 || this->data[k][i * c + j] > new_max )
                    throw(dimension_mismatch());
            }
        }
    }
}

Tensor Tensor::padding(int pad_h, int pad_w)const {
    int n_r = (r+2*pad_h), n_c = (c+2*pad_w);
    Tensor temp(n_r, n_c, d);

    for (int k = 0; k < d; k++)
    {
        for (int i = pad_h; i < (n_r-pad_h); i++)
        {
            for (int j = pad_w; j < (n_c-pad_w); j++)
            {
                temp.data[k][i * n_c + j] = data[k][(i-pad_h) * c + (j-pad_w)];
            }
        }
    }

    return temp;
}


Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end)const {
    if (col_start<0 || col_start>col_end || col_start>c || col_end<0 || col_end>c)
        throw(dimension_mismatch());
    else if (row_start<0 || row_start>col_end || row_start>r || row_end<0 || row_end>r)
        throw(dimension_mismatch());
    else if (depth_start<0 || depth_start>col_end || depth_start>c || depth_end<0 || depth_end>c)
        throw(dimension_mismatch());

    Tensor temp(row_end, col_end, depth_end);

    for (int k = depth_start; k < depth_end; k++){
        for (int i = row_start; i < row_end; i++){
            for (int j = col_start; j < col_end; j++){
                temp.data[k][i * col_end + j] = this->data[k][i * c + j];
            }
        }
    }

    temp.d = depth_end;
    temp.r = row_end;
    temp.c = col_end;

    return temp;
}

Tensor Tensor::concat(const Tensor &rhs, int axis=0)const {
    return *this;
}


Tensor Tensor::convolve(const Tensor &f)const {
    if (f.c % 2 == 0 || f.r % 2 == 0)
        throw (filter_odd_dimensions());
    if (f.d != this->d)
        throw (concat_wrong_dimension());
    
    Tensor new_padd_tensor;
    new_padd_tensor = this->padding( (f.r-1)/2 , (f.c-1)/2 );
    
    Tensor new_conv_tensor(this->r, this->c, this->d);

    for (int k = 0; k < new_conv_tensor.d; k++) {
        for (int i = 0; i <= new_conv_tensor.r; i++) {
            for (int j = 0; j <= new_conv_tensor.c; j++) {
                
                float tot = 0;
                for (int i_f=0; i_f<f.r; i++) {
                    for (int j_f=0; j_f<f.c; j++) {
                        tot = tot + ( new_padd_tensor.data[k][i_f*f.c+j_f] * f.data[k][i_f*f.c+j_f] );
                    }
                }
                new_conv_tensor.data[k][i*new_conv_tensor.c+j]=tot;

            }
        }
    }

    new_conv_tensor.clamp(0,255);

    new_conv_tensor.rescale();

    return new_conv_tensor;
}

int Tensor::rows()const {
    return r;
}

int Tensor::cols()const {
    return c;
}

int Tensor::depth()const {
    return d;
}
    


float Tensor::getMin(int k)const {
    float min = 0;
    if (k<d && k>=0) {
        min = data[k][0];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                if (data[k][i*c+j] < min)
                    min = data[k][i*c+j];
            }
        }
    } else {
        throw(index_out_of_bound());
    }
    return min;
}

float Tensor::getMax(int k)const {
    float max = 0;
    if (k<d && k>=0) {
        max = data[k][0];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                if (data[k][i*c+j] > max)
                    max = data[k][i*c+j];
            }
        }
    } else {
        throw(index_out_of_bound());
    }
    return max;
}


void Tensor::showSize()const {
    std::cout << r << " x " << c << " x " << d <<std::endl;
}
    
/*
friend ostream& Tensor::operator<< (ostream& stream, const Tensor & obj) {
    return *this;
}*/

 
void Tensor::read_file(string filename) {
    

}

 
void Tensor::write_file(string filename) {

}