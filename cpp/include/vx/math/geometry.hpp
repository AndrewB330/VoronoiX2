#ifndef VORONOIX_POLYGON_HPP
#define VORONOIX_POLYGON_HPP

#include "vec.hpp"
#include "mat.hpp"
#include <algorithm>
#include <optional>
#include <vector>

#define VX_ASSERT(exp, msg)
//if (!(exp)) std::cerr << msg << ":" << __FILE__ << ":" << __LINE__ << std::endl
// #define VX_LOG() std::cout << __FILE__ << ":" << __LINE__ << std::endl;

namespace vx {

    template<typename T>
    struct GeometryContext {
        T eps;
        T tight_eps;

        explicit GeometryContext(T maximal_span = 1.0) :
                eps(std::numeric_limits<T>::epsilon() * 1024 * 128 * maximal_span),
                tight_eps(std::numeric_limits<T>::epsilon() * 256 * maximal_span) {}
    };

    template<typename T, size_t DIM>
    struct HalfSpace {
        GeometryContext<T> ctx;
        Vec <T, DIM> normal;
        T offset = T(0);

        HalfSpace() = default;

        HalfSpace(const Vec <T, DIM> &normal, T offset, const GeometryContext<T> &ctx)
                : normal(normal), offset(offset), ctx(ctx) {}

        T orientedDistance(const Vec <T, DIM> &v) const {
            return dot(v, normal) - offset;
        }

        bool isAbove(const Vec <T, DIM> &v) const {
            return dot(v, normal) > offset + ctx.eps;
        }

        bool isBelow(const Vec <T, DIM> &v) const {
            return dot(v, normal) + ctx.eps < offset;
        }

        bool isOn(const Vec <T, DIM> &v) const {
            T dot_result = dot(v, normal);
            return dot_result <= offset + ctx.eps && dot_result + ctx.eps >= offset;
        }

        void flip() {
            normal = -normal;
            offset = -offset;
        }
    };

    template<typename T, size_t DIM>
    HalfSpace<T, DIM> makeHalfSpace(const std::vector<vx::Vec<T, DIM>> &points, const GeometryContext<T> &ctx) {
        VX_ASSERT(points.size() == DIM, "Num of points must be equal to DIM")
            << "Points: " << points.size() << " DIM: " << DIM << std::endl;
        Vec<T, DIM> c;
        for (auto v : points) c += v;
        c /= points.size();

        Mat<T, DIM, DIM> mat;
        for (size_t i = 0; i < DIM; i++) {
            mat.set_col(i, points[i] - c);
        }
        T offset = det(mat);
        Vec<T, DIM> normal;
        Mat<T, DIM, DIM> mat_copy = mat;
        for (size_t i = 0; i < DIM; i++) {
            std::fill_n(mat[i], DIM, T(1));
            normal[i] = det(mat);
            std::copy_n(mat_copy[i], DIM, mat[i]);
        }
        T d = length(normal);
        offset += dot(c, normal);
        normal /= d;
        offset /= d;
        for (const auto &p : points) {
            if (!(dot(normal, p) <= offset + ctx.tight_eps &&
                  dot(normal, p) + ctx.tight_eps >= offset)) {
                auto dd = dot(normal, p) - offset;
                std::cerr << "xx" << std::endl;
            }
            VX_ASSERT(
                    dot(normal, p) <= offset + ctx.tight_eps &&
                    dot(normal, p) + ctx.tight_eps >= offset,
                    "Floating point error is too big."
            ) << dot(normal, p) - offset << " eps: " << ctx.tight_eps << std::endl
              << mat << std::endl;
        }
        return HalfSpace(normal, offset, ctx);
    }

    template<size_t DIM>
    struct VerticesList {
        size_t vertices[DIM];
    };

    template<typename T, size_t DIM>
    struct Polyhedron {
        GeometryContext<T> ctx;
        std::vector<Vec < T, DIM>> points;
        std::vector<VerticesList<DIM>> facets;
    };


} // namespace vx

#endif //VORONOIX_POLYGON_HPP
