//  dtw.hpp
//  libDTW
//  Author: Gagandeep Singh, 10-02-2019
//

#ifndef dtw_hpp
#define dtw_hpp

#include <string>
#include <cmath>
#include <vector>
#include <iostream>

/**
 DTW path is returned using this struct
 */
typedef struct DTW_path {
    float* distances;
    int length;
} DTW_path;


/**
 Template for matrix
 */
template <class T>
class Matrix {
    int width = 0;
    int height = 0;
    T* data = NULL;
public:
    /**
     constructor for undefined initialization
     @param height number of rows
     @param width number of columns
     */
    Matrix(int height, int width) : height(height), width(width) {
        data = new T[width * height];
    }
    
    /**
     constructor to initialize all values to a value
     @param height number of rows
     @param width number of columns
     @param value value with which to initialize
     */
    Matrix(int height, int width, T value) : height(height), width(width) {
        data = new T[width * height];
        for (int i = 0; i < width * height; ++i)
            data[i] = value;
    }
    
    /**
     Destructor
     */
    ~Matrix(){
        delete[] data;
    }
    
    /**
     access elements by operator overload. Can be used to read and write
     @param r row
     @param c column
     */
    T& operator() (int r, int c){
        if (r >= height)
            throw std::out_of_range("row value out of bound");
        if (c >= width)
            throw std::out_of_range("column value out of bound");
        return data[r*width + c];
    }
    
    /**
     access const elements by operator overload, only used to read
     @param r row
     @param c column
     */
    T operator() (int r, int c) const{
        if (r >= height)
            throw std::out_of_range("row value out of bound");
        if (c >= width)
            throw std::out_of_range("column value out of bound");
        return data[r*width + c];
    }
    
    /**
     access number of columns
     */
    int get_num_cols() const{
        return width;
    }
    
    /**
     access number of rows
     */
    int get_num_rows() const{
        return height;
    }
    
    /**
     access data to read it
     */
    const T* Data() const {
        return data;
    }
};


/**
 Distance calculation function is implemented using this abstract class
 Currently only eucledean distance is implemented.
 */
class distance_func {
public:
    /**
     @param x array containing first vector
     @param y array containing second vector
     @param dim dimension(size) of both the vectors
     @return distance between these two vectors
     */
    virtual float operator() (const float* x, const float* y, int dim) = 0;
};


/**
 class for implementing euclidean distance between two vectors.
 */
class euclidean_dist : public distance_func {
public:
    /**
     eucledean distance between 2 vectors
     @param x 1st vector
     @param y 2nd vector
     @param dim dimension of each vector
     */
    virtual float operator() (const float* x, const float* y, int dim) {
        float d = 0;
        for (int i=0; i<dim; ++i)
            d += pow(x[i] - y[i], 2.0);
        return sqrt(d);
    }
};

/**
 Interface function for this C++ library. Returns a DTW_path struct given two
 sequences of vectors, the distance function, weight to put on diagonal movement
 and band width to use for the path.
 @param seq1 matrix filled with sequence 1
 @param seq2 matrix filled with sequence 2
 @param fn distance function
 @param w_diag weight on the diagonal
 @param band_win window to be considered along diagonal for DTW path
 @return dtw path
 */
DTW_path get_dtw_path(const Matrix<float>& seq1,
                      const Matrix<float>& seq2,
                      distance_func &fn, float w_diag, int band_win);

#endif /* dtw_hpp */

