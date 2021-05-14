#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "../dais_exc.h"

#define TYPE float

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

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

TYPE **padding(TYPE **&data, int d, int r, int c, int pad_h, int pad_w){
    int n_r = (r+2*pad_h), n_c = (c+2*pad_w);
    TYPE **temp = static_matrix_create(d, n_r, n_c);

    cout << n_r << " - " << n_c << endl << endl;

    for (int k = 0; k < d; k++)
    {
        for (int i = pad_h; i < (n_r-pad_h); i++)
        {
            for (int j = pad_w; j < (n_c-pad_w); j++)
            {
                temp[k][i * n_c + j] = data[k][(i-pad_h) * c + (j-pad_w)];
            }
        }
    }

    return temp;
}

void write_file(string filename, TYPE **&data, int d, int r, int c) {
    fstream f;
    f.open(filename, ios::out);
    int k, i, j;
    if(f.is_open()){
        f << d << endl << r << endl << c << endl;
        for (k = 0; k < d; k++){
            for (i = 0; i < r; i++){
                for (j = 0; j < c; j++){
                    f << data[k][i * c + j] << " ";
                }
                f << endl;
            }
        }
    }

    f.close();
}

class Tensor
{
private:

    float ** data = nullptr; //<-- you are free to change this data structure (don't use std::vectors or std::array)

    int r = 0;  // number of rows
    int c = 0;  // number of columns
    int d = 0;  // tensor depth

public:

    /**
     * Class constructor
     * 
     * Parameter-less class constructor 
     */
    Tensor();

    /**
     * Class constructor
     * 
     * Creates a new tensor of size r*c*d initialized at value v
     * 
     * @param r
     * @param c
     * @param d
     * @param v
     * @return new Tensor
     */
    Tensor(int r, int c, int d, float v = 0.0);

    /**
     * Class distructor
     * 
     * Cleanup the data when deallocated
     */
    ~Tensor();

    /**
     * Operator overloading ()
     * 
     * if indexes are out of bound throw index_out_of_bound() exception
     * 
     * @return the value at location [i][j][k]
     */
    float operator()(int i, int j, int k) const;

    /**
     * Operator overloading ()
     * 
     * Return the pointer to the location [i][j][k] such that the operator (i,j,k) can be used to 
     * modify tensor data.
     * 
     * If indexes are out of bound throw index_out_of_bound() exception
     * 
     * @return the pointer to the location [i][j][k]
     */
    float &operator()(int i, int j, int k);

    /**
     * Copy constructor
     * 
     * This constructor copies the data from another Tensor
     *      
     * @return the new Tensor
     */
    Tensor(const Tensor& that);

    /**
     * Operator overloading -
     * 
     * It performs the point-wise difference between two Tensors.
     * 
     * result(i,j,k)=this(i,j,k)-rhs(i,j,k)
     * 
     * The two tensors must have the same size otherwise throw a dimension_mismatch()
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator-(const Tensor &rhs)const;
    
     /**
     * Operator overloading +
     * 
     * It performs the point-wise sum between two Tensors.
     * 
     * result(i,j,k)=this(i,j,k)+rhs(i,j,k)
     * 
     * The two tensors must have the same size otherwise throw a dimension_mismatch()
     * 
     * @return returns a new Tensor containing the result of the operation
    */
    Tensor operator +(const Tensor &rhs)const;

    /**
     * Operator overloading *
     * 
     * It performs the point-wise product between two Tensors.
     * 
     * result(i,j,k)=this(i,j,k)*rhs(i,j,k)
     * 
     * The two tensors must have the same size otherwise throw a dimension_mismatch()
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator*(const Tensor &rhs)const;
    
    /**
     * Operator overloading /
     * 
     * It performs the point-wise division between two Tensors.
     * 
     * result(i,j,k)=this(i,j,k)/rhs(i,j,k)
     * 
     * The two tensors must have the same size otherwise throw a dimension_mismatch()
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator/(const Tensor &rhs)const;

    /**
     * Operator overloading - 
     * 
     * It performs the point-wise difference between a Tensor and a constant
     * 
     * result(i,j,k)=this(i,j,k)-rhs
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator-(const float &rhs)const;

    /**
     * Operator overloading +
     * 
     * It performs the point-wise sum between a Tensor and a constant
     * 
     * result(i,j,k)=this(i,j,k)+rhs
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator+(const float &rhs)const;

    /**
     * Operator overloading *
     * 
     * It performs the point-wise product between a Tensor and a constant
     * 
     * result(i,j,k)=this(i,j,k)*rhs
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator*(const float &rhs)const;

    /**
     * Operator overloading / between a Tensor and a constant
     * 
     * It performs the point-wise division between a Tensor and a constant
     * 
     * result(i,j,k)=this(i,j,k)/rhs
     * 
     * @return returns a new Tensor containing the result of the operation
     */
    Tensor operator/(const float &rhs)const;

    /**
     * Operator overloading = (assignment) 
     * 
     * Perform the assignment between this object and another
     * 
     * @return a reference to the receiver object
     */
    Tensor & operator=(const Tensor &other);

    /**
     * Random Initialization
     * 
     * Perform a random initialization of the tensor
     * 
     * @param mean The mean
     * @param std  Standard deviation
     */
    void init_random(float mean=0.0, float std=1.0);

    /**
     * Constant Initialization
     * 
     * Perform the initialization of the tensor to a value v
     * 
     * @param r The number of rows
     * @param c The number of columns
     * @param d The depth
     * @param v The initialization value
     */
    void init(int r, int c, int d, float v=0.0);

    /**
     * Tensor Clamp
     * 
     * Clamp the tensor such that the lower value becomes low and the higher one become high.
     * 
     * @param low Lower value
     * @param high Higher value 
     */
    void clamp(float low, float high);

    /**
     * Tensor Rescaling
     * 
     * Rescale the value of the tensor following this rule:
     * 
     * newvalue(i,j,k) = ((data(i,j,k)-min(k))/(max(k)-min(k)))*new_max
     * 
     * where max(k) and min(k) are the maximum and minimum value in the k-th channel.
     * 
     * new_max is the new value for the maximum
     * 
     * @param new_max New maximum vale
     */
    void rescale(float new_max=1.0);

    /**
     * Tensor padding
     * 
     * Zero pad a tensor in height and width, the new tensor will have the following dimensions:
     * 
     * (rows+2*pad_h) x (cols+2*pad_w) x (depth) 
     * 
     * @param pad_h the height padding
     * @param pad_w the width padding
     * @return the padded tensor
     */
    Tensor padding(int pad_h, int pad_w)const;

    /**
     * Subset a tensor
     * 
     * retuns a part of the tensor having the following indices:
     * row_start <= i < row_end  
     * col_start <= j < col_end 
     * depth_start <= k < depth_end
     * 
     * The right extrema is NOT included
     * 
     * @param row_start 
     * @param row_end 
     * @param col_start
     * @param col_end
     * @param depth_start
     * @param depth_end
     * @return the subset of the original tensor
     */
    Tensor subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end)const;

    /** 
     * Concatenate 
     * 
     * The function concatenates two tensors along a give axis
     * 
     * Example: this is of size 10x5x6 and rhs is of 25x5x6
     * 
     * if concat on axis 0 (row) the result will be a new Tensor of size 35x5x6
     * 
     * if concat on axis 1 (columns) the operation will fail because the number 
     * of rows are different (10 and 25).
     * 
     * In order to perform the concatenation is mandatory that all the dimensions 
     * different from the axis should be equal, other wise throw concat_wrong_dimension(). 
     *  
     * @param rhs The tensor to concatenate with
     * @param axis The axis along which perform the concatenation 
     * @return a new Tensor containing the result of the concatenation
     */
    Tensor concat(const Tensor &rhs, int axis=0)const;


    /** 
     * Convolution 
     * 
     * This function performs the convolution of the Tensor with a filter.
     * 
     * The filter f must have odd dimensions and same depth. 
     * 
     * Remeber to apply the padding before running the convolution
     *  
     * @param f The filter
     * @return a new Tensor containing the result of the convolution
     */
    Tensor convolve(const Tensor &f)const;

    /* UTILITY */

    /** 
     * Rows 
     * 
     * @return the number of rows in the tensor
     */
    int rows()const;

    /** 
     * Cols 
     * 
     * @return the number of columns in the tensor
     */
    int cols()const;

    /** 
     * Depth 
     * 
     * @return the depth of the tensor
     */
    int depth()const;
    
    /** 
     * Get minimum 
     * 
     * Compute the minimum value considering a particular index in the third dimension
     * 
     * @return the minimum of data( , , k)
     */
    float getMin(int k)const;

    /** 
     * Get maximum 
     * 
     * Compute the maximum value considering a particular index in the third dimension
     * 
     * @return the maximum of data( , , k)
     */
    float getMax(int k)const;

    /** 
     * showSize
     * 
     * shows the dimensions of the tensor on the standard output.
     * 
     * The format is the following:
     * rows" x "colums" x "depth
     * 
     */
    void showSize()const;
    
    /* IOSTREAM */

    /**
     * Operator overloading <<
     * 
     * Use the overaloading of << to show the content of the tensor.
     * 
     * You are free to chose the output format, btw we suggest you to show the tensor by layer.
     * 
     * [..., ..., 0]
     * [..., ..., 1]
     * ...
     * [..., ..., k]
     */
    friend ostream& operator<< (ostream& stream, const Tensor & obj);

    /**
     * Reading from file
     * 
     * Load the content of a tensor from a textual file.
     * 
     * The file should have this structure: the first three lines provide the dimensions while 
     * the following lines contains the actual data by channel.
     * 
     * For example, a tensor of size 4x3x2 will have the following structure:
     * 4
     * 3
     * 2
     * data(0,0,0)
     * data(0,1,0)
     * data(0,2,0)
     * data(1,0,0)
     * data(1,1,0)
     * .
     * .
     * .
     * data(3,1,1)
     * data(3,2,1)
     * 
     * if the file is not reachable throw unable_to_read_file()
     * 
     * @param filename the filename where the tensor is stored
     */
    void read_file(string filename);

    /**
     * Write the tensor to a file
     * 
     * Write the content of a tensor to a textual file.
     * 
     * The file should have this structure: the first three lines provide the dimensions while 
     * the following lines contains the actual data by channel.
     * 
     * For example, a tensor of size 4x3x2 will have the following structure:
     * 4
     * 3
     * 2
     * data(0,0,0)
     * data(0,1,0)
     * data(0,2,0)
     * data(1,0,0)
     * data(1,1,0)
     * .
     * .
     * .
     * data(3,1,1)
     * data(3,2,1)
     * 
     * @param filename the filename where the tensor should be stored
     */
    void write_file(string filename);

};

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

    if (k < 0 || k > d)
        throw(index_out_of_bound());
    else if (i < 0 || i > r)
        throw(index_out_of_bound());
    else if (j < 0 || j > c)
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
                temp.data[k][i * c + j] = (rhs.data[k][i * c + j]) / (this->data[k][i * c + j]);
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
    if (col_start>col_end || col_start>unsigned(c)  || col_end>unsigned(c) )
        throw(dimension_mismatch());
    else if (row_start>col_end || row_start>unsigned(r) || row_end>unsigned(r))
        throw(dimension_mismatch());
    else if (depth_start>col_end || depth_start>unsigned(c) || depth_end>unsigned(c) )
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

    else if(axis == 1 && this->c!=rhs.c)
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
    new_padd_tensor = this->padding( (f.r-1)/2 , (f.c-1)/2 );
    
    Tensor new_conv_tensor(this->r, this->c, this->d);

    for (int k = 0; k < new_conv_tensor.d; k++) {
        for (int i = 0; i < new_conv_tensor.r; i++) {
            for (int j = 0; j < new_conv_tensor.c; j++) {
                
                float tot = 0;
                for (int i_f=0; i_f<f.r; i++) {
                    for (int j_f=0; j_f<f.c; j++) {
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

        for (int k = 0; k < d; k++)
        {
            for (int i = 0; i < r; i++)
            {
                getline(f, q);
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
        f << d << "\n" << r << "\n" << c << "\n";
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



int main(){

    Tensor a(3, 3, 3, 5);
    Tensor b(a);
    cout << a << b;

    return 0;
}