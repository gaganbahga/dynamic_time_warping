//  python_wrapper.cpp
//  libDTW
//
//  Author: Gagandeep Singh, 11-02-2019
//

#include "dtw.hpp"

extern "C" {
    /**
     Get DTW path between two vectors as a struct
     @param seq1 first sequence of vectors
     @param seq2 second sequence of vectors
     @param len_seq1 length of first sequence
     @param len_seq2 length of second sequence
     @param dim dimension of the vectors in each sequence
     @param w_diag weight for diagonal movement in the dtw path
     @param band_win width of the band along diagonal where dtw path is allowed
     @return dtw path
     */
    DTW_path dtw_path(const float seq1[], const float seq2[], int len_seq1,
                      int len_seq2, int dim, float w_diag, int band_win){
        Matrix<float> m_seq1(len_seq1, dim);
        for (int i=0; i < len_seq1; ++i)
            for(int j=0; j < dim; ++j)
                m_seq1(i, j) = seq1[i * dim + j];
        
        Matrix<float> m_seq2(len_seq2, dim);
        for (int i=0; i < len_seq2; ++i)
            for(int j=0; j < dim; ++j)
                m_seq2(i, j) = seq2[i * dim + j];

        euclidean_dist ED;
        return get_dtw_path(m_seq1, m_seq2, ED, w_diag, band_win);
    }
    
    /**
     dtw_path function returns a struct which contains array. In order to deallocate
     memory of that array, call this function
     @param path pointer to the dtw path struct
     */
    void deallocate_DTW(DTW_path* path){
        delete[] path->distances;
    }
}
