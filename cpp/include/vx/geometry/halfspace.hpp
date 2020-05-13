#ifndef VORONOIX_HALFSPACE_HPP
#define VORONOIX_HALFSPACE_HPP

#include <vx/geometry/vector.hpp>
#include <vx/geometry/matrix.hpp>
#include <vx/utils.hpp>

namespace vx {

    template<typename T, size_t DIM>
    struct HalfSpace {
        vx::Vec<T, DIM> normal;
        T offset = T(0);

        HalfSpace() = default;

        HalfSpace(const Vec<T, DIM> &normal, T offset)
                : normal(normal), offset(offset) {}

        T orientedDistance(const Vec<T, DIM> &v) const {
            return dot(v, normal) - offset;
        }

        bool isAbove(const Vec<T, DIM> &v) const {
            return dot(v, normal) > offset + epsilon<T>();
        }

        bool isBelow(const Vec<T, DIM> &v) const {
            return dot(v, normal) + epsilon<T>() < offset;
        }

        bool isOn(const Vec<T, DIM> &v) const {
            T dot_result = dot(v, normal);
            return dot_result <= offset + epsilon<T>() && dot_result + epsilon<T>() >= offset;
        }

        void flip() {
            normal = -normal;
            offset = -offset;
        }
    };

    template<typename T, size_t DIM>
    vx::HalfSpace<T, DIM> makeHalfSpace(const std::vector<vx::Vec<T, DIM>> &points) {
        VX_ASSERT(points.size() == DIM, "Num of points must be equal to DIM")
            << "Points: " << points.size() << " DIM: " << DIM << std::endl;

        vx::Vec<T, DIM> c;
        for (const auto &p : points) c += p;
        c /= points.size();

        vx::Mat<T, DIM, DIM> mat;
        for (size_t i = 0; i < DIM; i++) {
            mat.set_col(i, points[i] - c);
        }

        T offset = det(mat);
        vx::Vec<T, DIM> normal;

        vx::Mat<T, DIM, DIM> mat_copy = mat;
        for (size_t i = 0; i < DIM; i++) {
            std::fill_n(mat[i], DIM, T(1));
            normal[i] = det(mat);
            std::copy_n(mat_copy[i], DIM, mat[i]);
        }

        offset += dot(c, normal);

        T d = length(normal);
        normal /= d;
        offset /= d;

        for (const auto &p : points) {
            VX_ASSERT(
                    dot(normal, p) <= offset + epsilon<T>() / 4 &&
                    dot(normal, p) + epsilon<T>() / 4 >= offset,
                    "Floating point error is too big."
            ) << dot(normal, p) - offset << " eps: " << epsilon<T>() / 4 << std::endl
              << mat << std::endl;
        }
        return vx::HalfSpace<T, DIM>(normal, offset);
    }

    template<typename T, size_t DIM, typename TList>
    vx::HalfSpace<T, DIM> makeHalfSpace(const std::vector<vx::Vec<T, DIM>> &points, const TList &vertices) {
        std::vector<vx::Vec<T, DIM>> selected;
        selected.reserve(vertices.size());
        for (size_t i : vertices) selected.push_back(points[i]);
        return makeHalfSpace(selected);
    }

}

#endif //VORONOIX_HALFSPACE_HPP
