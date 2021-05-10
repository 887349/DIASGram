#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"
#include "vector_lib\static_array.h"

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
    data = static_matrix_create(0,0,0,0.0);

    r = 0; // number of rows
    c = 0; // number of columns
    d = 0; // tensor depth
}

Tensor::Tensor(int r, int c, int d, float v=0.0) {

}

Tensor::Tensor(const Tensor& that) {

}

Tensor::~Tensor() {
    static_del(data, d);
}

float Tensor::operator()(int i, int j, int k) const {
    return 0;
}

float& Tensor::operator()(int i, int j, int k) {
    float* a;
    return *a;
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
    
/*
    K È DEEP (D)
    PER ACCEDERE A UNA CELLA NELLA DEEP K FACCIO :
        data[i*c+j][k] o data[k][i*c+j]
    SECONDO ME QUELLO CHE HO SCRITTO NEL CODICE MA NON SONO SICURO 
*/
/*float Tensor::getMin(int k)const {
    float min = 0;data[k][0];
    if (k<d) {
        min = data[k][0];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                if (data[k][i*c+j] < min)
                    min = data[k][i*c+j];
            }
        }
    }
    return min;
}*/

/*float Tensor::getMax(int k)const {
    float max = 0;
    if (k<d) {
        max = data[k][0];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                if (data[k][i*c+j] > max)
                    max = data[k][i*c+j];
            }
        }
    }
    return max;
}*/

void Tensor::showSize()const {
    std::cout << r << " x " << c << " x " << d <<std::endl;
}
    
/*friend ostream& Tensor::operator<< (ostream& stream, const Tensor & obj) {
    return *this;
}*/

 
void Tensor::read_file(string filename) {

}

 
void Tensor::write_file(string filename) {

}