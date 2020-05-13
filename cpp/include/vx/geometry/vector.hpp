#ifndef VORONOIX_VECTOR_HPP
#define VORONOIX_VECTOR_HPP

#include <iostream>
#include <optional>
#include <cmath>

namespace vx {

    template<typename T, size_t DIM>
    struct Vec {

        T data[DIM] = {};

        Vec() {}

        const T &operator[](size_t i) const {
            return data[i];
        }

        T &operator[](size_t i) {
            return data[i];
        }
    };

    template<typename T>
    struct Vec<T, 2> {
        T x, y;

        explicit Vec(T x = T(), T y = T()) : x(x), y(y) {}

        const T &operator[](size_t i) const {
            return (i ? y : x);
        }

        T &operator[](size_t i) {
            return (i ? y : x);
        }
    };

    template<typename T>
    struct Vec<T, 3> {
        T x, y, z;

        explicit Vec(T x = T(), T y = T(), T z = T()) : x(x), y(y), z(z) {}

        const T &operator[](size_t i) const {
            return (i ? (i > 1 ? z : y) : x);
        }

        T &operator[](size_t i) {
            return (i ? (i > 1 ? z : y) : x);
        }
    };

    template<typename T, size_t DIM>
    Vec<T, DIM> &operator+=(Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        for (size_t i = 0; i < DIM; i++) {
            a[i] += b[i];
        }
        return a;
    }

    template<typename T, size_t DIM>
    Vec<T, DIM> operator+(Vec<T, DIM> a, const Vec<T, DIM> &b) {
        return a += b;
    }

    template<typename T, size_t DIM>
    Vec<T, DIM> &operator-=(Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        for (size_t i = 0; i < DIM; i++) {
            a[i] -= b[i];
        }
        return a;
    }

    template<typename T, size_t DIM>
    Vec<T, DIM> operator-(Vec<T, DIM> a, const Vec<T, DIM> &b) {
        return a -= b;
    }

    template<typename T, size_t DIM, typename K>
    Vec<T, DIM> &operator*=(Vec<T, DIM> &a, K k) {
        for (size_t i = 0; i < DIM; i++) {
            a[i] *= k;
        }
        return a;
    }

    template<typename T, size_t DIM, typename K>
    Vec<T, DIM> operator*(Vec<T, DIM> a, K k) {
        return a *= k;
    }

    template<typename T, size_t DIM, typename K>
    Vec<T, DIM> &operator/=(Vec<T, DIM> &a, K k) {
        for (size_t i = 0; i < DIM; i++) {
            a[i] /= k;
        }
        return a;
    }

    template<typename T, size_t DIM, typename K>
    Vec<T, DIM> operator/(Vec<T, DIM> a, K k) {
        return a /= k;
    }

    template<typename T, size_t DIM>
    Vec<T, DIM> operator-(const Vec<T, DIM> &a) {
        return Vec<T, DIM>() - a;
    }

    template<typename T, size_t DIM>
    bool operator<(const Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        for (size_t i = 0; i < DIM; i++) {
            if (a[i] != b[i]) return a[i] < b[i];
        }
        return false;
    }

    template<typename T, size_t DIM>
    bool operator==(const Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        for (size_t i = 0; i < DIM; i++) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    template<typename T, size_t DIM>
    T length_square(const Vec<T, DIM> &v) {
        T res = T();
        for (size_t i = 0; i < DIM; i++) res += v[i] * v[i];
        return res;
    }

    template<typename T, size_t DIM>
    T length(const Vec<T, DIM> &v) {
        return std::sqrt(length_square(v));
    }

    template<typename T, size_t DIM>
    T dist_square(const Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        return length_square(a - b);
    }

    template<typename T, size_t DIM>
    T dist(const Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        return length(a - b);
    }

    template<typename T, size_t DIM>
    T dot(const Vec<T, DIM> &a, const Vec<T, DIM> &b) {
        T res = T();
        for (size_t i = 0; i < DIM; i++) res += a[i] * b[i];
        return res;
    }

    template<typename T, size_t DIM>
    std::ostream &operator<<(std::ostream &out, const Vec<T, DIM> &v) {
        out << '(';
        for (size_t i = 0; i < DIM; i++)
            out << v[i] << (i + 1 < DIM ? ", " : ")");
        return out;
    }

    template<typename T>
    inline T cross_product(const vx::Vec<T, 2> &a, const vx::Vec<T, 2> &b) {
        return a.x * b.y - a.y * b.x;
    }

    template<typename T>
    inline T cross_product(const vx::Vec<T, 2> &a, const vx::Vec<T, 2> &b, const vx::Vec<T, 2> &c) {
        return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
    }

    template<typename T>
    bool is_inside(const vx::Vec<T, 2> &a, const vx::Vec<T, 2> &b, const vx::Vec<T, 2> &c, const vx::Vec<T, 2> &o) {
        return cross_product(b - a, o - a) > 0 &&
               cross_product(c - b, o - b) > 0 &&
               cross_product(a - c, o - c) > 0;
        // TODO: check epsilon
    }

    template<typename T>
    vx::Vec<T, 3> cross_product(const vx::Vec<T, 3> &a, const vx::Vec<T, 3> &b) {
        return vx::Vec<T, 3>(a.y * b.z - a.z * b.y, a.x * b.z - a.z * b.x, a.x * b.y - a.y * b.x);
    }

    template<typename T, size_t DIM>
    struct CmpReversed {
        bool operator()(const vx::Vec<T, DIM> &a, const vx::Vec<T, DIM> &b) const {
            for (size_t i = DIM; i > 0; i--) {
                if (a[i - 1] != b[i - 1]) return a[i - 1] < b[i - 1];
            }
            return false;
        }
    };

} // namespace vx

#endif //VORONOIX_VECTOR_HPP
