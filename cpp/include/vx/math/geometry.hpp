#ifndef VORONOIX_POLYGON_HPP
#define VORONOIX_POLYGON_HPP

#include <vx/math/vec.hpp>
#include <optional>
#include <vector>

namespace vx {

    template<typename T, size_t DIM>
    struct Simplex {
        Vec<T, DIM + 1> points;

        Simplex() {}
    };

    template<typename T, size_t DIM>
    struct Plane {
        Vec<T, DIM> points;

        Plane() {}
    };

    template<typename T, size_t DIM>
    struct Ridge {

    };


} // namespace vx

#endif //VORONOIX_POLYGON_HPP
