#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"
#include "matrix_lib\matrix.h"

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
                data[k][i * 2 + j] = that.data[k][i * 2 + j];
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

    if(i < d || i < d)
        throw(index_out_of_bound());
    else if (j < r || j < r)
        throw(index_out_of_bound());
    else if (i < c || i < c)
        throw(index_out_of_bound());
    else{
        
        res = data[k][(i*c)+j];
    }

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
    return rhs;
}

Tensor Tensor::operator +(const Tensor &rhs)const {
    return rhs;
}

Tensor Tensor::operator*(const Tensor &rhs)const {
    return rhs;
}

Tensor Tensor::operator/(const Tensor &rhs)const {
    return rhs;
}

Tensor Tensor::operator-(const float &rhs)const {
    return *this;
}

Tensor Tensor::operator+(const float &rhs)const {
    return *this;
}

Tensor Tensor::operator*(const float &rhs)const {
    return *this;
}

Tensor Tensor::operator/(const float &rhs)const {
    return *this;
}

Tensor & Tensor::operator=(const Tensor &other) {
    return *this;
}

void Tensor::init_random(float mean=0.0, float std=1.0) {

}

void Tensor::init(int r, int c, int d, float v=0.0) {

}

void Tensor::clamp(float low, float high) {

}

void Tensor::rescale(float new_max=1.0) {

}

Tensor Tensor::padding(int pad_h, int pad_w)const {
    return *this;
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end)const {
    return *this;
}

Tensor Tensor::concat(const Tensor &rhs, int axis=0)const {
    return *this;
}


Tensor Tensor::convolve(const Tensor &f)const {
    return f;
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
    

/* ADD THROW ERROR
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
*/

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