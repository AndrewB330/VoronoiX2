#ifndef VORONOIX_ALG_CONVEX_HULL_HPP
#define VORONOIX_ALG_CONVEX_HULL_HPP

#include "alg_convex_hull_2d.hpp"
#include "alg_convex_hull_kd.hpp"

namespace vx {

    template<typename T, size_t DIM>
    class ConvexHull {
    public:
        explicit ConvexHull(const std::vector<vx::Vec<T, DIM>> &points);

        std::vector<std::array<size_t, DIM>> getFacets();

        std::vector<vx::HalfSpace<T, DIM>> getHalfSpaces();

        size_t getSize();

        bool verify(const std::vector<vx::Vec<T, DIM>> &points);

    private:
        vx::QuickHull<T, DIM> quick_hull;
    };

    template<typename T>
    class ConvexHull<T, 2> {
    public:
        explicit ConvexHull(const std::vector<vx::Vec<T, 2>> &points);

        std::vector<std::array<size_t, 2>> getFacets();

        std::vector<vx::HalfSpace<T, 2>> getHalfSpaces();

        size_t getSize();

        bool verify(const std::vector<vx::Vec<T, 2>> &points);

    private:
        vx::ConvexHull2D<T> convex_hull;
        std::vector<Vec<T, 2>> points;
    };

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    ConvexHull<T, DIM>::ConvexHull(const std::vector<vx::Vec<T, DIM>> &points):quick_hull(points) {}

    template<typename T, size_t DIM>
    std::vector<std::array<size_t, DIM>> ConvexHull<T, DIM>::getFacets() {
        return quick_hull.getFacets();
    }

    template<typename T, size_t DIM>
    bool ConvexHull<T, DIM>::verify(const std::vector<Vec<T, DIM>> &points) {
        return quick_hull.isBounded() && quick_hull.isConvex();
    }

    template<typename T, size_t DIM>
    size_t ConvexHull<T, DIM>::getSize() {
        return quick_hull.getFacets().size(); // TODO: optimize
    }

    template<typename T, size_t DIM>
    std::vector<HalfSpace<T, DIM>> ConvexHull<T, DIM>::getHalfSpaces() {
        return quick_hull.getHalfSpaces();
    }

    template<typename T>
    ConvexHull<T, 2>::ConvexHull(const std::vector<Vec<T, 2>> &points):convex_hull(GrahamScan(points)), points(points) {}


    template<typename T>
    std::vector<std::array<size_t, 2>> ConvexHull<T, 2>::getFacets() {
        return convex_hull.getSegments();
    }

    template<typename T>
    std::vector<HalfSpace<T, 2>> ConvexHull<T, 2>::getHalfSpaces() {
        std::vector<HalfSpace<T, 2>> hs;
        for(auto & facet : convex_hull.getSegments()) {
            auto dif = points[facet[1]] - points[facet[0]];
            auto len = length(dif);
            auto normal = Vec<T, 2>(-dif.y, dif.x) / len;
            hs.emplace_back(normal, len);
        }
        return hs;
    }

    template<typename T>
    bool ConvexHull<T, 2>::verify(const std::vector<vx::Vec<T, 2>> &points) {
        return convex_hull.verifyResult(points);
    }

    template<typename T>
    size_t ConvexHull<T, 2>::getSize() {
        return convex_hull.getSize();
    }

} // namespace vx

#endif //VORONOIX_ALG_CONVEX_HULL_HPP
