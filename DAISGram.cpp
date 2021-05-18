#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram(){

}

DAISGram::~DAISGram(){
    
}

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){ 
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);    
            data(i,j,2) = (float) img.blue_at(j,i);   
        }                
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 * 
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));                   
        }                
    }

    img.write(filename);

}


/**
 * Generate Random Image
 * 
 * Generate a random image from nois
 * 
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */  
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}

int DAISGram::getRows(){
    return data.rows();
}

int DAISGram::getCols(){
    return data.cols();
}

int DAISGram::getDepth(){
    return data.depth();
}

DAISGram DAISGram::brighten(float bright){
    DAISGram res;
    res.data = data;
    res.data = res.data + bright;
    res.data.clamp(0, 255);

    return res;
}

DAISGram DAISGram::grayscale(){
    if(getDepth() <= 0)
        throw(index_out_of_bound());

    DAISGram res;
    res.data = data;
    
    float depth = float(getDepth());
        
    for (int i=0; i<getRows(); i++){
        for(int j=0; j<getCols(); j++){
            float somma = 0;
            for(int k=0; k<getDepth(); k++) {
                somma+=data(i, j, k);
            }
            float media = somma/depth;
            for(int k=0; k<depth; k++) {
                res.data(i, j, k)=media;
            }
        }
    }

    return res;
}

DAISGram DAISGram::warhol(){
    DAISGram res;
    Tensor tdx, bsx, bdx; 
    tdx = data;
    bsx = data;
    bdx = data;

    tdx.swap_channel(0, 1);
    bsx.swap_channel(2, 1);
    bdx.swap_channel(0, 2);

    res.data = (data.concat(tdx, 1)).concat((bsx.concat(bdx, 1)), 0);

    return res;
}

DAISGram DAISGram::sharpen(){
    
    Tensor filter(3,3,data.depth());
    
    for(int k=0; k<filter.depth(); k++) {

        filter(0,0,k)=filter(0,2,k)=filter(2,0,k)=filter(2,2,k)=0;
        filter(0,1,k)=filter(1,0,k)=filter(1,2,k)=filter(2,1,k)=-1;
        filter(1,1,k)=5;
    }

    DAISGram res;
    res.data = data.convolve(filter);
    
    res.data.clamp(0, 255); //NON DOVREBBE SERVIRE PERCHE VIENE GIA EFFETTUATO IN CONVOLVE
    
    return res;
}

DAISGram DAISGram::emboss(){
    Tensor filter(3,3,this->data.depth());

    for(int k=0; k<filter.depth(); k++) {
        
        filter(0,0,k) = -2;
        filter(0,1,k) = filter(1,0,k) = -1;
        filter(2,0,k) = filter(0,2,k) = 0;
        filter(1,1,k) = filter(2,1,k) = filter(0,2,k) = 1;
        filter(2,2,k) = 2;
    }
    
    DAISGram res;
    res.data = this->data.convolve(filter);
    
    res.data.clamp(0, 255);
    
    return res;
}

DAISGram DAISGram::smooth(int h){
    DAISGram res;
    if(h <= 0)
        throw(index_out_of_bound());
    
    float a = (float)1/(h*h);

    Tensor filter(h, h, this->data.depth(), a);

    res.data = this->data.convolve(filter);

    return res;
}

DAISGram DAISGram::edge(){
    
    Tensor filter(3, 3, this->data.depth(), -1);
    
    for(int k=1; k<filter.depth(); k++) {
        filter(1,1,k)=8;
    }

    DAISGram res = this->grayscale();


    res.data = res.data.convolve(filter);

    res.data.clamp(0, 255); //NON DOVREBBE SERVIRE PERCHE VIENE GIA EFFETTUATO IN CONVOLVE
    
    return res;
}


DAISGram DAISGram::blend(const DAISGram & rhs, float alpha){

    if (this->getRows()!=rhs->getRows() or this->getRows()!=rhs->getRows())
        throw(dimension_mismatch);

    DAISGram res;

    res.data = alpha*this.data + (1-alpha)*rhs.data;
    // Non sono sicura se basti cosi o se serva un ciclo
    return res;
}


DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){
    return *this;
}

DAISGram DAISGram::equalize(){
    return *this;
}
