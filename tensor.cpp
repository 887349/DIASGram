#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "tensor.h"
#include "matrix.cpp"

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

Tensor::Tensor(int r, int c, int d, float v) {
    data = static_matrix_create(d, r, c);
    matrix_init(data, d, r, c, v);

    this->r = r; // number of rows
    this->c = c; // number of columns
    this->d = d; // tensor depth
}

Tensor::Tensor(const Tensor& that) {

    data = static_matrix_create(that.d, that.r, that.c);
    /* *this = that; */

    r = that.r;
    c = that.c;
    d = that.d;

    for (int k = 0; k < d; k++)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                this->data[k][i * c + j] = (that.data[k][i * c + j]);
            }
        }
    }   
}

Tensor::~Tensor() {
    matrix_del(data, d);
}

float Tensor::operator()(int i, int j, int k) const {

    if(k < d || k < d)
        throw(index_out_of_bound());
    if (i < r || i < r)
        throw(index_out_of_bound());
    if (j < c || j < c)
        throw(index_out_of_bound());
    else
        return data[k][(i*c)+j];
}

float& Tensor::operator()(int i, int j, int k) {
    
    if (k < 0 || k > d)
        throw(index_out_of_bound());
    if (i < 0 || i > r)
        throw(index_out_of_bound());
    if (j < 0 || j > c)
        throw(index_out_of_bound());
    else
        return (data[k][(i * c) + j]);
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
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) - (rhs.data[k][i * c + j]);
            }
        }
    }

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
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) / (rhs.data[k][i * c + j]);
            }
        }
    }

    return temp;
}

Tensor Tensor::operator-(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) - rhs;
            }
        }
    }

    return temp;
}

Tensor Tensor::operator+(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) + rhs;
            }
        }
    }

    return temp;
}

Tensor Tensor::operator*(const float &rhs)const {

    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) * rhs;
            }
        }
    }

    return temp;
}

Tensor Tensor::operator/(const float &rhs)const {
    Tensor temp(r, c, d);

    for (int k = 0; k < d; k++){
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temp.data[k][i * c + j] = (this->data[k][i * c + j]) / rhs;
            }
        }
    }

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
                this->data[k][i * c + j] = (other.data[k][i * c + j]);
            }
        }
    }
    
    return *this;
}

void Tensor::init(int r, int c, int d, float v) {

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

void Tensor::rescale(float new_max) {
    for (int k = 0; k < d; k++)
    {
        float min = getMin(k);
        float max = getMax(k);
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                if (max==min)
                    data[k][i * c + j] = new_max;
                else
                    data[k][i * c + j] = ((data[k][i * c + j] - min) / (max - min)) * new_max;

                
                if (data[k][i * c + j] < 0 || data[k][i * c + j] > new_max )
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
    if (col_start>col_end || col_start>unsigned(c)  || col_end>unsigned(c) )
        throw(dimension_mismatch());
    if (row_start>col_end || row_start>unsigned(r) || row_end>unsigned(r))
        throw(dimension_mismatch());
    if (depth_start>col_end || depth_start>unsigned(d) || depth_end>unsigned(d) )
        throw(dimension_mismatch());


    Tensor temp(row_end, col_end, depth_end);

    for (unsigned int k = depth_start; k < depth_end; k++){
        for (unsigned int i = row_start; i < row_end; i++){
            for (unsigned int j = col_start; j < col_end; j++){
                temp.data[k][i * col_end + j] = this->data[k][i * c + j];
            }
        }
    }

    return temp;
}

Tensor Tensor::concat(const Tensor &rhs, int axis)const {
    if(axis == 0 && this->r!=rhs.r)
        throw concat_wrong_dimension();

    if(axis == 1 && this->c!=rhs.c)
        throw concat_wrong_dimension();

    if(axis == 0){
        int new_r = this->r + rhs.r;
        Tensor temp(new_r, c, d);

        for (int k = 0; k < d; k++)
        {
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                     temp.data[k][i * c + j] = data[k][i * c + j];
                }
            }
        }

        for (int k = 0; k < d; k++)
        {
            for (int i = r; i < new_r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    int app_i=i-r;
                    temp.data[k][i * c + j] = rhs.data[k][app_i * c + j];
                }
            }
        }
        return temp;

    }
    else if (axis==1) {
        int new_c = this->c + rhs.c;
        Tensor temp(r, new_c, d);
        for (int k = 0; k < d; k++){
            for (int i = 0; i < temp.r; i++){
                for (int j = 0; j < temp.c; j++){
                    if (j<this->c) 
                        temp.data[k][i * temp.c + j] = this->data[k][i * this->c + j];
                    else {
                        int app_j = j-this->c;
                        temp.data[k][i * temp.c + j] = rhs.data[k][i * rhs.c + app_j];
                    }
                }
            }
        }
        return temp;
    }
    else {
        Tensor temp;
        return temp;
    }
}

Tensor Tensor::convolve(const Tensor &f)const {
    if (f.c % 2 == 0 || f.r % 2 == 0)
        throw (filter_odd_dimensions());
    if (f.d != this->d)
        throw (concat_wrong_dimension());
    
    Tensor new_padd_tensor;
    new_padd_tensor = this->padding( (f.r-1)/2 , (f.c-1)/2);


    Tensor new_conv_tensor(this->r, this->c, this->d);
    for (int k = 0; k < new_conv_tensor.d; k++) {
        for (int i = 0; i < new_conv_tensor.r; i++) {
            for (int j = 0; j < new_conv_tensor.c; j++) {
                
                float tot = 0;
                for (int i_f=0; i_f<f.r; i_f++) {
                    for (int j_f=0; j_f<f.c; j_f++) {
                        tot = tot + ( new_padd_tensor.data[k][(i+i_f)*new_padd_tensor.c+(j+j_f)] * f.data[k][i_f*f.c+j_f] );
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
    
ostream& operator<< (ostream& stream, const Tensor & obj){

    int k, i, j;

    stream << obj.d << "\n" << obj.r << "\n" << obj.c << "\n\n";
    for (k = 0; k < obj.d; k++)
    {
        for (i = 0; i < obj.r; i++)
        {
            for (j = 0; j < obj.c; j++)
            {
                stream << obj.data[k][i * obj.c + j] << " ";
            }
            stream << "\n";
        }
        stream << "\n";
    }

    return stream;
}

void Tensor::read_file(string filename) {
    ifstream f{filename};
    if(f.is_open()) {
        string q;
        getline(f, q);
        r = stoi(q);
        getline(f, q);
        c = stoi(q);
        getline(f, q);
        d = stoi(q);
        
        matrix_del(data, d);
        data = static_matrix_create ( d, r, c);

        for (int k = 0; k < d; k++)
        {
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    f >> data[k][i * c + j];
                }
            }
        }
    }    
    else
        throw(unable_to_read_file());

    f.close();
}

void Tensor::write_file(string filename) {
    fstream f;
    f.open(filename, ios::out);

    int k, i, j;
    if(f.is_open()){
        f << r << "\n" << c << "\n" << d << "\n";
        for (k = 0; k < d; k++)
        {
            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    f << data[k][i * c + j] << " ";
                }
                f << "\n";
            }
            f << "\n";
        }
    }
    else
        throw(unable_to_write_file());

    f.close();
}

void Tensor::swap_channel(int a, int b){
    float *temp;
    temp = data[a];
    data[a] = data[b];
    data[b] = temp;
}