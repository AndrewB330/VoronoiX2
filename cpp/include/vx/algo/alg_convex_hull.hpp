#ifndef VORONOIX_ALG_CONVEX_HULL_HPP
#define VORONOIX_ALG_CONVEX_HULL_HPP

#include "alg_convex_hull_2d.hpp"
#include "alg_convex_hull_kd.hpp"

namespace vx {

    template<typename T, size_t DIM>
    class ConvexHull {
    public:
        explicit ConvexHull(const std::vector<vx::Vec<T, DIM>> &points);

        std::vector<std::array<size_t, 2>> getFacets();

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

        size_t getSize();

        bool verify(const std::vector<vx::Vec<T, 2>> &points);

    private:
        vx::ConvexHull2D<T> convex_hull;
    };

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    ConvexHull<T, DIM>::ConvexHull(const std::vector<vx::Vec<T, DIM>> &points):quick_hull(points) {}

    template<typename T, size_t DIM>
    std::vector<std::array<size_t, 2>> ConvexHull<T, DIM>::getFacets() {
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

    template<typename T>
    ConvexHull<T, 2>::ConvexHull(const std::vector<Vec<T, 2>> &points):convex_hull(GrahamScan(points)) {}


    template<typename T>
    std::vector<std::array<size_t, 2>> ConvexHull<T, 2>::getFacets() {
        return convex_hull.getSegments();
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
