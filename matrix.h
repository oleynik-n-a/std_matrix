#ifndef MATRIX_H
#define MATRIX_H

#pragma once

#include <iostream>

namespace std {

class matrix_out_of_range : public std::out_of_range {
public:
    matrix_out_of_range() : std::out_of_range("Index out of range") {
    }
};

class matrix_not_square : public std::runtime_error {
public:
    matrix_not_square() : std::runtime_error("Matrix is not square") {
    }
};

class matrix_unexist : public std::runtime_error {
public:
    matrix_unexist() : std::runtime_error("Matrix unexist") {
    }
};

template <typename _Tp, size_t _M, size_t _N>
struct matrix {
    using value_type = _Tp;
    using size_type = size_t;
    using pointer = _Tp*;
    using const_pointer = const _Tp*;
    using reference = _Tp&;
    using const_reference = const _Tp&;
    using iterator = _Tp*;
    using const_iterator = const _Tp*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    value_type __elems_[_M][_N];

    size_t rows() const {
        return _M;
    }

    size_t columns() const {
        return _N;
    }

    value_type& at(size_t a, size_t b) {
        if (a >= _M || b >= _N) {
            throw matrix_out_of_range();
        }
        return __elems_[a][b];
    }

    const value_type& at(size_t a, size_t b) const {
        if (a >= _M || b >= _N) {
            throw matrix_out_of_range();
        }
        return __elems_[a][b];
    }

    value_type& operator()(size_t a, size_t b) {
        return __elems_[a][b];
    }
    
    const value_type& operator()(size_t a, size_t b) const {
        return __elems_[a][b];
    }

    value_type& front() {
        return &__elems_[0][0];
    }

    const value_type& front() const {
        return &__elems_[0][0];
    }

    value_type& back() {
        return &__elems_[_M - 1][_N - 1];
    }

    const value_type& back() const {
        return &__elems_[_M - 1][_N - 1];
    }

    value_type& data() {
        return &__elems_;
    }

    const value_type& data() const {
        return &__elems_;
    }

    iterator begin() noexcept {  // NOLINT
        return &__elems_[0][0];
    }

    iterator end() noexcept {  // NOLINT
        return &__elems_[_M - 1][_N];
    }

    const_iterator cbegin() const noexcept {  // NOLINT
        return &__elems_[0][0];
    }

    const_iterator cend() const noexcept {  // NOLINT
        return &__elems_[_M - 1][_N];
    }

    reverse_iterator rbegin() noexcept {  // NOLINT
        return std::reverse_iterator<iterator>(end());
    }

    reverse_iterator rend() noexcept {  // NOLINT
        return std::reverse_iterator<iterator>(begin());
    }

    const_reverse_iterator crbegin() const noexcept {  // NOLINT
        return std::reverse_iterator<iterator>(cend());
    }

    const_reverse_iterator crend() const noexcept {  // NOLINT
        return std::reverse_iterator<iterator>(cbegin());
    }

    bool empty() const {
        return _M == 0 || _N == 0;
    }

    size_t size() const {
        return _M * _N;
    }

    void fill(value_type value) {
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                __elems_[i][j] = value;
            }
        }
    }

    void swap(matrix<_Tp, _M, _N>& other) {
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                std::swap(__elems_[i][j], other.__elems_[i][j]);
            }
        }
    }

    value_type trace() const {
        value_type trace = __elems_[0][0];
        for (size_t i = 1; i < _M && i < _N; ++i) {
            trace += __elems_[i][i];
        }

        return trace;
    }
 
    value_type determinant() const {
        if (_M != _N) {
            throw matrix_not_square();
        }

        matrix<value_type, _N, _N> temp = *this;
        for (size_t i = 0; i < _N - 1; i++) {
            for (size_t j = i + 1; j < _N; j++) {
                for (size_t k = i + 1; k < _N; k++) {
                    temp(j, k) = (temp(j, k) * temp(i, i) - temp(j, i) * temp(i, k));
                    if (i) {
                        temp(j, k) /= temp(i - 1, i - 1);
                    }
                }
            }
        }

        return temp(_N - 1, _N - 1);
    }

    value_type minor(size_t m, size_t n) const {
        if (_M != _N) {
            throw matrix_not_square();
        }

        if (m >= _N || n >= _N || _N == 1) {
            throw std::matrix_out_of_range();
        }

        matrix<value_type, _N - 1, _N - 1> minor_matrix;
        bool is_i = false;
        for (size_t i = 0; i < _N - 1; ++i) {
            bool is_j = false;
            if (i == m) {
                is_i = true;
            }
            for (size_t j = 0; j < _N - 1; ++j) {
                if (j == n) {
                    is_j = true;
                }
                minor_matrix(i, j) = __elems_[is_i ? i + 1 : i][is_j ? j + 1 : j];
            }
        }

        for (size_t i = 0; i < _N - 2; i++) {
            for (size_t j = i + 1; j < _N - 1; j++) {
                for (size_t k = i + 1; k < _N - 1; k++) {
                    minor_matrix(j, k) = (minor_matrix(j, k) * minor_matrix(i, i) -
                                          minor_matrix(j, i) * minor_matrix(i, k));
                    if (i) {
                        minor_matrix(j, k) /= minor_matrix(i - 1, i - 1);
                    }
                }
            }
        }

        return minor_matrix(_N - 2, _N - 2);
    }

    matrix<value_type, _N, _M> get_transposed() const {
        matrix<value_type, _N, _M> result;
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                result(j, i) = __elems_[i][j];
            }
        }
        return result;
    }

    matrix<value_type, _M, _N> get_inversed() const {
        if (_M != _N) {
            throw matrix_not_square();
        }

        value_type det = determinant();
        if (det == 0) {
            throw matrix_unexist();
        }

        matrix<value_type, _M, _N> result;
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                result(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * minor(i, j);
            }
        }

        result /= det;

        return result.get_transposed();
    }

    matrix<value_type, _M, _N> operator+(const matrix<value_type, _M, _N>& other) const {
        matrix<value_type, _M, _N> result;
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                result.__elems_[i][j] = __elems_[i][j] + other.__elems_[i][j];
            }
        }
        return result;
    }

    matrix<value_type, _M, _N> operator-(const matrix<value_type, _M, _N>& other) const {
        matrix<value_type, _M, _N> result;
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                result.__elems_[i][j] = __elems_[i][j] - other.__elems_[i][j];
            }
        }
        return result;
    }

    template <size_t _K>
    matrix<value_type, _M, _K> operator*(const matrix<value_type, _N, _K>& other) const {
        matrix<value_type, _M, _K> result;
        result.fill(0);
        for (size_t i = 0; i < _M; ++i) {
            for (size_t k = 0; k < _K; ++k) {
                for (size_t j = 0; j < _N; ++j) {
                    result.__elems_[i][k] += __elems_[i][j] * other.__elems_[j][k];
                }
            }
        }
        return result;
    }

    matrix<value_type, _M, _N>& operator+=(const matrix<value_type, _M, _N>& other) {
        *this = *this + other;
        return *this;
    }

    matrix<value_type, _M, _N>& operator-=(const matrix<value_type, _M, _N>& other) {
        *this = *this - other;
        return *this;
    }

    template <size_t _K>
    matrix<value_type, _M, _K>& operator*=(const matrix<value_type, _N, _K>& other) {
        *this = *this * other;
        return *this;
    }

    bool operator==(const matrix<value_type, _M, _N>& other) const {
        for (size_t i = 0; i < _M; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                if (__elems_[i][j] != other.__elems_[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const matrix<value_type, _M, _N>& other) const {
        return !(*this == other);
    }
};

template <typename value_type, size_t _M, size_t _N>
matrix<value_type, _M, _N> operator*(const matrix<value_type, _M, _N>& matr, const value_type& value) {
    matrix<value_type, _M, _N> result;
    for (size_t i = 0; i < _M; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            result(i, j) = matr(i, j) * value;
        }
    }
    return result;
}

template <typename value_type, size_t _M, size_t _N>
matrix<value_type, _M, _N> operator*(const value_type& value, const matrix<value_type, _M, _N>& matr) {
    return matr * value;
}

template <typename value_type, size_t _M, size_t _N>
matrix<value_type, _M, _N> operator/(const matrix<value_type, _M, _N>& matr, const value_type& value) {
    matrix<value_type, _M, _N> result;
    for (size_t i = 0; i < _M; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            result(i, j) = matr(i, j) / value;
        }
    }
    return result;
}

template <typename value_type, size_t _M, size_t _N>
matrix<value_type, _M, _N>& operator*=(matrix<value_type, _M, _N>& matr, const value_type& value) {
    matr = matr * value;
    return matr;
}

template <typename value_type, size_t _M, size_t _N>
matrix<value_type, _M, _N>& operator/=(matrix<value_type, _M, _N>& matr, const value_type& value) {
    matr = matr / value;
    return matr;
}

template <typename value_type, size_t _M, size_t _N>
std::ostream& operator<<(std::ostream& ostream, matrix<value_type, _M, _N> matr) {
    for (size_t i = 0; i < _M; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            ostream << matr(i, j);
            if (j != _N - 1) {
                ostream << ' ';
            }
        }
        ostream << '\n';
    }

    return ostream;
}

template <typename value_type, size_t _M, size_t _N>
std::istream& operator>>(std::istream& istream, matrix<value_type, _M, _N>& matr) {
    for (size_t i = 0; i < _M; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            istream >> matr(i, j);
        }
    }

    return istream;
}
 
}  // namespace std

#endif
