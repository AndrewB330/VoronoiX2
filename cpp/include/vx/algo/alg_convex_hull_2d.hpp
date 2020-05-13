#ifndef VORONOIX_ALG_CONVEX_HULL_2D_HPP
#define VORONOIX_ALG_CONVEX_HULL_2D_HPP

#include <vx/geometry/vector.hpp>
#include <vector>
#include <array>

namespace vx {

    template<typename T>
    class ConvexHull2D {
    public:
        bool verifyResult(const std::vector<vx::Vec<T, 2>> &points);

        std::vector<size_t> getVertices();

        std::vector<std::array<size_t, 2>> getSegments();

        size_t getSize();

    protected:
        ConvexHull2D() = default;

        std::vector<size_t> result; // in counter-clockwise order
    };

    template<typename T>
    class GrahamScan : public ConvexHull2D<T> {
    public:
        explicit GrahamScan(const std::vector<vx::Vec<T, 2>> &points);
    };

    template<typename T>
    class GiftWrapping : public ConvexHull2D<T> {
    public:
        explicit GiftWrapping(const std::vector<vx::Vec<T, 2>> &points);
    };

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T>
    GrahamScan<T>::GrahamScan(const std::vector<vx::Vec<T, 2>> &points):ConvexHull2D<T>() {
        // TODO: verify input
        auto &result = this->result;
        std::vector<size_t> order;
        for (size_t i = 0; i < points.size(); i++) order.emplace_back(i);
        auto leftmost_index = std::min_element(points.begin(), points.end()) - points.begin();
        std::swap(order[0], order[leftmost_index]);
        auto pivot = points[order[0]];
        std::sort(order.begin() + 1, order.end(), [&](const auto &i, const auto &j) {
            T cross_results = cross_product(points[i] - pivot, points[j] - pivot);
            if (std::abs(cross_results) < vx::epsilon<T>()) {
                return dist_square(points[i], pivot) < dist_square(points[j], pivot);
            }
            return cross_results > 0;
        });

        result.push_back(order[0]);
        result.push_back(order[1]);

        for (size_t i = 2; i < order.size(); i++) {
            while (result.size() >= 2) {
                auto prev = points[result[result.size() - 1]];
                auto prev_prev = points[result[result.size() - 2]];
                if (cross_product(prev - prev_prev, points[order[i]] - prev_prev) <= 0) {
                    result.pop_back();
                } else {
                    break;
                }
            }
            result.push_back(order[i]);
        }

        VX_ASSERT(this->verifyResult(points), "Convex hull verification failed");

    }

    template<typename T>
    GiftWrapping<T>::GiftWrapping(const std::vector<vx::Vec<T, 2>> &points):ConvexHull2D<T>() {
        auto &result = this->result;
        auto current_index = std::min_element(points.begin(), points.end()) - points.begin();
        auto start_index = current_index;
        std::vector<char> used(points.size(), 0);
        do {
            size_t best_index = current_index;
            for (size_t i = 0; i < points.size(); i++) {
                if (!used[i] && i != current_index) {
                    if (best_index == current_index) {
                        best_index = i;
                    } else {
                        auto v1 = points[i] - points[current_index];
                        auto v2 = points[best_index] - points[current_index];
                        T cross_res = cross_product(v1, v2);
                        if (cross_res > 0) {
                            best_index = i;
                        } else if (cross_res + vx::epsilon<T>() >= 0) {
                            if (length(v1) > length(v2)) {
                                best_index = i;
                            }
                        }
                    }
                }
            }
            if (best_index == current_index) {
                std::cout << best_index << std::endl;
            }
            VX_ASSERT(best_index != current_index, "Not found");
            current_index = best_index;
            used[current_index] = true;
            result.push_back(current_index);
        } while (current_index != start_index);
        VX_ASSERT(this->verifyResult(points), "Convex hull verification failed");
    }

    template<typename T>
    bool ConvexHull2D<T>::verifyResult(const std::vector<vx::Vec<T, 2>> &points) {
        for (const auto &segment : getSegments()) {
            for (const auto &point : points) {
                auto v1 = points[segment[1]] - points[segment[0]];
                auto v2 = point - points[segment[0]];
                if (cross_product(v1, v2) + vx::epsilon<T>() < 0)
                    return false;
            }
        }
        return true;
    }

    template<typename T>
    std::vector<size_t> ConvexHull2D<T>::getVertices() {
        return result;
    }

    template<typename T>
    std::vector<std::array<size_t, 2>> ConvexHull2D<T>::getSegments() {
        std::vector<std::array<size_t, 2>> segments;
        for (size_t i = 0; i + 1 < result.size(); i++) {
            segments.push_back(std::array<size_t, 2>{result[i], result[i + 1]});
        }
        segments.push_back(std::array<size_t, 2>{result[result.size() - 1], result[0]});
        return segments;
    }

    template<typename T>
    size_t ConvexHull2D<T>::getSize() {
        return result.size();
    }

} // namespace vx

#endif //VORONOIX_ALG_CONVEX_HULL_2D_HPP
