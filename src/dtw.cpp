//  dtw.cpp
//  libDTW
//  Author: Gagandeep Singh, 11-02-2019
//

#include "dtw.hpp"

#define TOP 0
#define SIDE 1
#define DIAGONAL 2
#define NOT_VISITED -1

/**
 Distance between two sequences of vectors
 @param f1 pointer to the first feature array
 @param f2 pointer to the second feature array
 @param dim dimension of each feature vector
 @param pdist pointer to array where array is stored, memory is already allocated
 */
void pair_distance(const float* f1, const float* f2, int dim, int band_win, Matrix<float>& pdist, distance_func& d) {
    int f1_len = pdist.get_num_rows(),
    f2_len = pdist.get_num_cols(),
    y_min = 0,
    y_max = f2_len;
    for (int x = 0; x < f1_len; ++x){
        if (band_win > 1){
            y_min = std::ceil((float)x * f2_len/f1_len - band_win);
            y_max = std::ceil((float)x * f2_len/f1_len + band_win);
            y_min = y_min > -1 ? y_min : 0;
            y_max = y_max < f2_len ? y_max : f2_len;
        }
        for (int y = y_min; y < y_max; ++y){
            pdist(x, y) = d(f1 + x * dim, f2 + y * dim, dim);
        }
    }
}

/**
 Calculates min and argmin of 3 numbers, leaving any neg number. Results are passed by reference
 @param x 1st num
 @param y 2nd num
 @param z 3rd num
 @param _min minimum
 @param argmin position of min
 */
void min(float x, float y, float z, float &_min, int &argmin){
    if (x < 0.)
        x = std::numeric_limits<float>::max();
    if (y < 0.)
        y = std::numeric_limits<float>::max();
    if (z < 0.)
        z = std::numeric_limits<float>::max();
    if (x < y && x < z){
        _min = x;
        argmin = 0;
    } else if (y < x && y < z) {
        _min = y;
        argmin = 1;
    } else {
        _min = z;
        argmin = 2;
    }
}

void backtrack(const Matrix<float>& dist, const Matrix<char> &direction,
               std::vector<float>& values, std::vector<std::tuple<int, int>>& indices){
    int path_len = 0;
    std::vector<std::tuple<int, int>> reverse_path;
    reverse_path.reserve(dist.get_num_rows() + dist.get_num_cols());

    int x = dist.get_num_rows() - 1,
    y = dist.get_num_cols() - 1;
    reverse_path.push_back(std::make_tuple(x, y));
    ++path_len;
    while(x > 0 && y > 0){
        int d = round(direction(x, y));
        switch (d) {
            case TOP:
                --x;
                break;
            case SIDE:
                --y;
                break;
            case DIAGONAL:
                --x;
                --y;
                break;
            case NOT_VISITED:
                ;
        }
        reverse_path.push_back(std::make_tuple(x, y));
        ++path_len;
    }
    values.resize(path_len);
    indices.reserve(path_len);

    for (int i = 0; i < path_len; ++i){
        indices[path_len - i - 1] = reverse_path[i];
        values[path_len - i - 1] = dist(std::get<0>(reverse_path[i]), std::get<1>(reverse_path[i]));
    }
}

/**
 Computes DTW given distances matrix. The outcome is the DTW path values and path indices
 @param dist distance matrix
 @param path_val DTW path values
 @param path_indices DTW path indices
 @param w_diag weight put of diagonal steps
 @return total dtw cost
 */
float calculate_dtw(Matrix<float>& dist, std::vector<float>& path_val, std::vector<std::tuple<int, int>>& path_indices, float w_diag) {
    Matrix<float> alpha(dist.get_num_rows(), dist.get_num_cols());
    Matrix<char> direction(dist.get_num_rows(), dist.get_num_cols());
    alpha(0, 0) = dist(0, 0);

    // y == 0
    for (int i = 1; i < alpha.get_num_rows(); ++i){
        if (dist(i, 0) < 0.){ // if -1
            direction(i, 0) = NOT_VISITED;
            alpha(i, 0) = (float)NOT_VISITED;
            continue;
        }
        alpha(i, 0) = alpha(i - 1, 0) + dist(i, 0);
        direction(i, 0) = (float)TOP;
    }
    // x == 0
    for (int i = 1; i < alpha.get_num_cols(); ++i){
        if (dist(0, i) < 0.){
            direction(0, i) = (float)NOT_VISITED;
            alpha(0, i) = (float)NOT_VISITED;
            continue;
        }
        alpha(0, i) = alpha(0, i - 1) + dist(0, i);
        direction(0, i) = (float)SIDE;
    }
    // interior points
    for (int x = 1; x < alpha.get_num_rows(); ++x) {
        for (int y = 1; y < alpha.get_num_cols(); ++y) {
            if (dist(x, y) < 0.){
                direction(x, y) = (float)NOT_VISITED;
                alpha(x, y) = (float)NOT_VISITED;
                continue;
            }
            float s1 = alpha(x - 1, y), //top
            s2 = alpha(x, y - 1), // side
            s3 = alpha(x - 1, y - 1); // diagonal
            // alpha[x * cols + y] = (float) smin(s1, s2, s3, eta) + pdist[x * cols + y];
            float _min;
            int argmin;
            min(s1, s2, s3, _min, argmin);
            if (argmin == DIAGONAL)
                alpha(x, y) = _min + w_diag * dist(x, y);
            else
                alpha(x, y) = _min + dist(x, y);
            switch (argmin){
                case TOP:
                    direction(x, y) = (float)TOP;
                    break;
                case SIDE:
                    direction(x, y) = (float)SIDE;
                    break;
                case DIAGONAL:
                    direction(x, y) = (float)DIAGONAL;
            }
        }
    }
    backtrack(dist, direction, path_val, path_indices);

    return alpha(alpha.get_num_rows() - 1, alpha.get_num_cols() - 1);
}

DTW_path get_dtw_path(const Matrix<float>& seq1,
                      const Matrix<float>& seq2,
                      distance_func &fn, float w_diag, int band_win){
    assert(seq1.get_num_cols() == seq2.get_num_cols());
    Matrix<float> distances(seq1.get_num_rows(), seq2.get_num_rows(), -1);
    pair_distance(seq1.Data(), seq2.Data(), seq1.get_num_cols(), band_win, distances, fn);
    std::vector<float> path_values;
    std::vector<std::tuple<int, int>> path_indices;
    float total_cost = calculate_dtw(distances, path_values, path_indices, w_diag);
    float *data = new float[path_values.size()];
    for (int i = 0; i < path_values.size(); i++)
        data[i] = path_values[i];
    DTW_path path = {.distances = data, .length = (int)path_values.size()};
    return path;
}

