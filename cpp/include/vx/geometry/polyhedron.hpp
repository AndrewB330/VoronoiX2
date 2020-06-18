#ifndef VORONOIX_POLYGON_HPP
#define VORONOIX_POLYGON_HPP

#include <algorithm>
#include <optional>
#include <vector>

#include <vx/geometry/halfspace.hpp>
#include <vx/algo/alg_convex_hull.hpp>

#include "vector.hpp"
#include "matrix.hpp"

namespace vx {

    template<typename T, size_t DIM>
    struct Polyhedron {
        explicit Polyhedron() {}

        Polyhedron(const std::vector<Vec<T, DIM>> &points, const std::vector<std::array<size_t, DIM>> &facets);

        std::vector<Vec<T, DIM>> points;
        std::vector<std::array<size_t, DIM>> facets;
    };

    template<typename T, size_t DIM>
    vx::Polyhedron<T, DIM> makePolyhedron(const std::vector<vx::Vec<T, DIM>> &points);

    template<typename T, size_t DIM>
    vx::Polyhedron<T, DIM> makePolyhedron(const std::vector<vx::HalfSpace<T, DIM>> &hs, const vx::Vec<T, DIM> &inside);

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    vx::Polyhedron<T, DIM>::Polyhedron(const std::vector<Vec<T, DIM>> &points,
                                       const std::vector<std::array<size_t, DIM>> &facets)
            :points(points), facets(facets) {}

    template<typename T, size_t DIM>
    vx::Polyhedron<T, DIM> makePolyhedron(const std::vector<vx::Vec<T, DIM>> &points) {
        vx::ConvexHull<T, DIM> hull(points);
        return vx::Polyhedron<T, DIM>(points, hull.getFacets());
    }

    template<typename T, size_t DIM>
    vx::Polyhedron<T, DIM>
    makePolyhedron(const std::vector<vx::HalfSpace<T, DIM>> &hs_, const vx::Vec<T, DIM> &inside) {
        // point-plane duality
        // Point = Normal / Offset

        std::vector<vx::Vec<T, DIM>> dual_vertices;
        for (const auto &hs: hs_) {
            VX_ASSERT(hs.isBelow(inside), "Given point isn't below all halfspaces");
            auto offset = hs.offset - dot(inside, hs.normal);
            dual_vertices.push_back(hs.normal / offset);
        }
        vx::ConvexHull<T, DIM> dual_hull(dual_vertices);

        std::vector<Vec<T, DIM>> result_vertices;
        for (const auto &hs: dual_hull.getHalfSpaces()) {
            result_vertices.push_back(hs.normal / hs.offset + inside);
        }

        return makePolyhedron(result_vertices);
    }

} // namespace vx

#endif //VORONOIX_POLYGON_HPP
