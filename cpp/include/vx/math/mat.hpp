#ifndef VORONOIX_MAT_HPP
#define VORONOIX_MAT_HPP

#include <iostream>
#include <cassert>
#include <memory>

#include "vec.hpp"

namespace vx {

    template<typename T, size_t ROWS, size_t COLS>
    struct Mat {
        T data[ROWS][COLS] = {};

        Mat() {}

        void load_identity() {
            assert(ROWS == COLS);
            for (size_t i = 0; i < ROWS; i++) {
                for (size_t j = 0; j < COLS; j++) {
                    data[i][j] = (i == j ? T(1) : T());
                }
            }
        }

        void transpose() {
            assert(ROWS == COLS);
            for (size_t i = 0; i < ROWS; i++) {
                for (size_t j = 0; j < i; j++) {
                    std::swap(data[i][j], data[j][i]);
                }
            }
        }

        const T *operator[](size_t i) const {
            return data[i];
        }

        T *operator[](size_t i) {
            return data[i];
        }

        Vec <T, COLS> get_row(size_t i) const {
            Vec<T, COLS> row;
            for (size_t j = 0; j < COLS; j++) row[j] = data[i][j];
            return row;
        }

        Vec <T, ROWS> get_col(size_t j) const {
            Vec<T, ROWS> col;
            for (size_t i = 0; i < ROWS; i++) col[i] = data[i][j];
            return col;
        }

        void set_row(size_t i, const Vec <T, COLS> &row) {
            for (size_t j = 0; j < COLS; j++) data[i][j] = row[j];
        }

        void set_col(size_t j, const Vec <T, ROWS> &col) {
            for (size_t i = 0; i < ROWS; i++) data[i][j] = col[i];
        }

        void swap_rows(size_t i1, size_t i2) {
            for (size_t j = 0; j < COLS; j++) {
                std::swap(data[i1][j], data[i2][j]);
            }
        }
    };

    template<typename T, size_t DIM1, size_t DIM2, size_t DIM3>
    Mat<T, DIM1, DIM3> operator*(const Mat<T, DIM1, DIM2> &a, const Mat<T, DIM2, DIM3> &b) {
        Mat<T, DIM1, DIM3> r;
        for (int i = 0; i < DIM1; i++) {
            for (int k = 0; k < DIM2; k++) {
                for (int j = 0; j < DIM3; j++) {
                    r[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return r;
    }

    template<typename T, size_t R, size_t C>
    std::ostream &operator<<(std::ostream &out, const Mat<T, R, C> &mat) {
        for (size_t i = 0; i < R; i++) {
            out << mat.get_row(i) << std::endl;
        }
        return out;
    }

    template<typename T, size_t DIM>
    T det(const Mat<T, DIM, DIM> &mat) {
        T eps = 1e-08;
        T res = T(1);
        auto det_mat = mat;
        for (size_t i = 0; i < DIM; ++i) {
            size_t pivot = i;
            T abs_max = std::abs(det_mat[i][i]);
            for (size_t j = i + 1; j < DIM; j++) {
                T abs_ji = std::abs(det_mat[j][i]);
                if (abs_max < abs_ji) {
                    abs_max = abs_ji;
                    pivot = j;
                }
            }
            if (!(eps < abs_max)) {
                return T(0);
            }
            if (pivot != i) {
                res = -res;
                det_mat.swap_rows(i, pivot);
            }
            T diag_val = det_mat[i][i];
            res *= diag_val; // det is multiple of diagonal elements
            for (size_t j = i + 1; j < DIM; j++) {
                T &mat_ji = det_mat[j][i];
                mat_ji /= diag_val;
                for (size_t k = i + 1; k < DIM; k++) {
                    det_mat[j][k] -= mat_ji * det_mat[i][k];
                }
            }
        }
        return res;
    }

    template<typename T>
    T det(const Mat<T, 2, 2> &mat) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    template<typename T>
    T det(const Mat<T, 3, 3> &mat) {
        return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) +
               mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) +
               mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    }

}

#endif //VORONOIX_MAT_HPP
