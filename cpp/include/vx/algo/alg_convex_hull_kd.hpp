#ifndef VORONOIX_ALG_CONVEX_HULL_KD_HPP
#define VORONOIX_ALG_CONVEX_HULL_KD_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <ctime>

#include <vx/geometry/polyhedron.hpp>
#include <vx/geometry/halfspace.hpp>

namespace vx {

    // TODO: remove shared_pointer to increase performance
    template<typename T, size_t DIM>
    struct Facet {
        using FacetPtr = std::shared_ptr<Facet<T, DIM>>;
        std::array<size_t, DIM> vertices = {};
        FacetPtr neighbours[DIM] = {};

        size_t furthest = 0;
        T furthest_dist = T(0);

        std::vector<size_t> outside;
        //std::vector<size_t> coplanar;

        vx::HalfSpace<T, DIM> half_space = {};

        explicit Facet() = default;

        void replaceNeighbour(FacetPtr old, FacetPtr current);

        void updateHalfSpace(const std::vector<vx::Vec<T, DIM>> &points);

        void checkAndPutVertex(size_t i, const Vec<T, DIM> &v);

        void flip();

        T dist(const Vec<T, DIM> &v) const;

        bool isOutside(const Vec<T, DIM> &v) const;

        bool isCoplanar(const Vec<T, DIM> &v) const;

        bool operator<(const Facet<T, DIM> &b);
    };

    template<typename T, size_t DIM>
    struct Ridge {
        std::array<size_t, DIM - 1> vertices;

        Ridge(const std::shared_ptr<Facet<T, DIM>> &facet, size_t side);

        bool operator==(const Ridge &r) const;
    };


    template<typename T, size_t DIM>
    struct RidgeHasher {
        size_t operator()(const Ridge<T, DIM> &r) const;
    };

    template<typename T, size_t DIM>
    class QuickHull {
    private:
        vx::Vec<T, DIM> inside;
        std::vector<vx::Vec<T, DIM>> points;
        std::set<std::shared_ptr<vx::Facet<T, DIM>>> facet_candidates;
        std::unordered_set<std::shared_ptr<vx::Facet<T, DIM>>> facets;

        void addFacet(const std::shared_ptr<vx::Facet<T, DIM>> &facet);

        void deleteFacet(const std::shared_ptr<vx::Facet<T, DIM>> &facet);

        void connectAllNewFacets(const std::vector<std::shared_ptr<vx::Facet<T, DIM>>> &new_facets);

        vx::VMat<T> buildGramianMatrix(const std::vector<size_t> &vertices, size_t other);

        size_t getFurthestVertex(const std::vector<size_t> &vertices);

        void buildInitialFacets();

        void buildConvexHull();

    public:
        explicit QuickHull(const std::vector<vx::Vec<T, DIM>> &points);

        std::vector<vx::HalfSpace<T, DIM>> getHalfSpaces();

        std::vector<std::array<size_t, DIM>> getFacets();

        bool isConvex() const;

        bool isBounded() const;
    };

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    void vx::Facet<T, DIM>::replaceNeighbour(vx::Facet<T, DIM>::FacetPtr old, vx::Facet<T, DIM>::FacetPtr current) {
        auto it = std::find(neighbours, neighbours + DIM, old);
        VX_ASSERT(it != neighbours + DIM, "Can't find facet in neighbour_facets");
        *it = current;
    }

    template<typename T, size_t DIM>
    void vx::Facet<T, DIM>::updateHalfSpace(const std::vector<Vec<T, DIM>> &points) {
        half_space = makeHalfSpace(points, vertices);
    }

    template<typename T, size_t DIM>
    void vx::Facet<T, DIM>::checkAndPutVertex(size_t i, const vx::Vec<T, DIM> &v) {
        /*if (isCoplanar(v)) { // TODO:
            coplanar.push_back(i);
        } else */
        if (isOutside(v)) {
            outside.push_back(i);
            T cur_dist = dist(v);
            if (cur_dist > furthest_dist) {
                furthest_dist = cur_dist;
                furthest = i;
            }
        }
    }

    template<typename T, size_t DIM>
    void vx::Facet<T, DIM>::flip() {
        std::swap(vertices[0], vertices[1]);
        std::swap(neighbours[0], neighbours[1]);
        half_space.flip();
    }

    template<typename T, size_t DIM>
    T vx::Facet<T, DIM>::dist(const vx::Vec<T, DIM> &v) const {
        return half_space.orientedDistance(v);
    }

    template<typename T, size_t DIM>
    bool vx::Facet<T, DIM>::isOutside(const vx::Vec<T, DIM> &v) const {
        return half_space.isAbove(v);
    }

    template<typename T, size_t DIM>
    bool vx::Facet<T, DIM>::isCoplanar(const vx::Vec<T, DIM> &v) const {
        return half_space.isOn(v);
    }

    template<typename T, size_t DIM>
    vx::Ridge<T, DIM>::Ridge(const std::shared_ptr<Facet<T, DIM>> &facet, size_t side) {
        for (size_t i = 0; i < DIM; i++) {
            if (i != side) {
                vertices[i - (i > side)] = facet->vertices[i];
            }
        }
        std::sort(vertices.begin(), vertices.end());
    }

    template<typename T, size_t DIM>
    bool vx::Ridge<T, DIM>::operator==(const vx::Ridge<T, DIM> &r) const {
        return std::equal(vertices.begin(), vertices.end(), r.vertices.begin());
    }


    template<typename T, size_t DIM>
    size_t vx::RidgeHasher<T, DIM>::operator()(const vx::Ridge<T, DIM> &r) const {
        size_t h = 0;
        for (size_t i = 0; i < DIM - 1; i++) {
            h ^= (r.vertices[i] * 1583) ^ (r.vertices[i]);
        }
        return h;
    }

    template<typename T, size_t DIM>
    bool vx::Facet<T, DIM>::operator<(const vx::Facet<T, DIM> &other) {
        return furthest_dist > other.furthest_dist;
    }

    template<typename T, size_t DIM>
    void vx::QuickHull<T, DIM>::addFacet(const std::shared_ptr<Facet<T, DIM>> &facet) {
        facets.insert(facet);
        if (!facet->outside.empty())
            facet_candidates.insert(facet);
    }

    template<typename T, size_t DIM>
    void vx::QuickHull<T, DIM>::deleteFacet(const std::shared_ptr<Facet<T, DIM>> &facet) {
        facet_candidates.erase(facet);
        facets.erase(facet);
        for (size_t i = 0; i < DIM; i++) {
            facet->neighbours[i] = nullptr;
        }
    }

    template<typename T, size_t DIM>
    void vx::QuickHull<T, DIM>::connectAllNewFacets(const std::vector<std::shared_ptr<Facet<T, DIM>>> &new_facets) {
        using RidgeInfo = std::pair<std::shared_ptr<Facet<T, DIM>>, size_t>;
        std::unordered_map<vx::Ridge<T, DIM>, RidgeInfo, vx::RidgeHasher<T, DIM>> ridge_map{};
        for (const auto &facet : new_facets) {
            for (size_t side = 0; side < DIM; side++) {
                auto ridge = vx::Ridge<T, DIM>(facet, side);
                if (ridge_map.find(ridge) != ridge_map.end()) {
                    auto[other_facet, other_side] = ridge_map[ridge];
                    other_facet->neighbours[other_side] = facet;
                    facet->neighbours[side] = other_facet;
                } else {
                    ridge_map[ridge] = {facet, side};
                }
            }
        }
    }

    template<typename T, size_t DIM>
    vx::VMat<T> vx::QuickHull<T, DIM>::buildGramianMatrix(const std::vector<size_t> &vertices, size_t other) {
        VMat<T> gramian_matrix(vertices.size(), vertices.size());
        for (size_t row = 0; row < gramian_matrix.rows; row++) {
            for (size_t col = 0; col < gramian_matrix.cols; col++) {
                gramian_matrix[row][col] = dot(
                        points[vertices[row]] - points[other],
                        points[vertices[col]] - points[other]
                );
            }
        }
        return gramian_matrix;
    }

    template<typename T, size_t DIM>
    size_t vx::QuickHull<T, DIM>::getFurthestVertex(const std::vector<size_t> &vertices) {
        T max_volume = T(0);
        size_t index = -1;
        for (size_t i = 0; i < points.size(); i++) {
            if (std::find(vertices.begin(), vertices.end(), i) != vertices.end()) continue;
            T volume = T(0);
            if (vertices.empty()) {
                for (size_t dim = 0; dim < DIM; dim++)
                    volume = std::max(volume, std::abs(points[i][dim]));
            } else {
                volume = det<T>(buildGramianMatrix(vertices, i));
            }
            if (volume > max_volume) {
                max_volume = volume;
                index = i;
            }
        }
        VX_ASSERT(index < points.size(), "Out of range") << "Index: " << index << std::endl;
        return index;
    }

    template<typename T, size_t DIM>
    void vx::QuickHull<T, DIM>::buildInitialFacets() {
        std::vector<size_t> simplex;
        for (size_t i = 0; i < DIM + 1; i++) simplex.push_back(getFurthestVertex(simplex));

        for (size_t v : simplex) inside += points[v];
        inside /= simplex.size();

        std::vector<std::shared_ptr<Facet<T, DIM>>> simplex_facets(DIM + 1);
        for (size_t i = 0; i < DIM + 1; i++) {
            simplex_facets[i] = std::make_shared<Facet<T, DIM>>();
        }
        for (size_t k = 0; k < DIM + 1; k++) {
            auto facet = simplex_facets[k];
            for (size_t i = 0; i < DIM + 1; i++) {
                if (i != k) {
                    facet->vertices[i - (i > k)] = simplex[i];
                    facet->neighbours[i - (i > k)] = simplex_facets[i];
                }
            }
            facet->updateHalfSpace(points);
            if (facet->isOutside(inside)) {
                facet->flip();
            }
            for (size_t i = 0; i < points.size(); i++)
                facet->checkAndPutVertex(i, points[i]);
            addFacet(facet);
        }
    }

    template<typename T, size_t DIM>
    void vx::QuickHull<T, DIM>::buildConvexHull() {
        std::unordered_set<size_t> outside_vertex_candidates{};
        while (!facet_candidates.empty()) {
            auto facet = *facet_candidates.begin();
            outside_vertex_candidates.clear();
            size_t furthest_vertex = facet->furthest;

            std::vector<std::shared_ptr<Facet<T, DIM>>> visible_set;

            std::unordered_map<std::shared_ptr<Facet<T, DIM>>, bool> visited = {};
            std::queue<std::shared_ptr<Facet<T, DIM>>> queue;

            visited[facet] = true;
            queue.push(facet);

            std::vector<std::pair<std::shared_ptr<Facet<T, DIM>>, size_t>> ridges;

            // BFS on visible facets
            while (!queue.empty()) {
                auto current_facet = queue.front();
                queue.pop();
                visible_set.push_back(current_facet);
                for (auto v : current_facet->outside) outside_vertex_candidates.insert(v);
                for (size_t i = 0; i < DIM; i++) {
                    auto neighbour = current_facet->neighbours[i];
                    if (!neighbour->isOutside(points[furthest_vertex])) {
                        ridges.emplace_back(current_facet, i);
                    } else if (!visited[neighbour]) {
                        queue.push(neighbour);
                        visited[neighbour] = true;
                    }
                }
            }

            std::vector<std::shared_ptr<Facet<T, DIM>>> new_facets;

            for (auto[facet, connection] : ridges) {
                auto old_facet = facet;
                auto new_facet = std::make_shared<Facet<T, DIM>>();

                std::copy_n(old_facet->vertices.begin(), DIM, new_facet->vertices.begin());
                new_facet->vertices[connection] = furthest_vertex;
                new_facet->neighbours[connection] = old_facet->neighbours[connection];
                new_facet->neighbours[connection]->replaceNeighbour(old_facet, new_facet);

                new_facet->updateHalfSpace(points);

                new_facet->outside.reserve(outside_vertex_candidates.size() / 2);
                for (size_t i : outside_vertex_candidates)
                    new_facet->checkAndPutVertex(i, points[i]);
                // TODO: should we check only candidates?not all?

                new_facets.push_back(new_facet);
                addFacet(new_facet);
                VX_ASSERT(!new_facet->isOutside(inside), "New facet should be directed outside");
            }

            connectAllNewFacets(new_facets);

            for (const auto &to_delete: visible_set) deleteFacet(to_delete);
        }
    }

    template<typename T, size_t DIM>
    vx::QuickHull<T, DIM>::QuickHull(const std::vector<vx::Vec<T, DIM>> &points) : points(points) {
        VX_ASSERT(points.size() > DIM, "Number of points should be at least DIM+1 for vx::QuickHull algrithm");

        buildInitialFacets();
        buildConvexHull();

        VX_ASSERT(isBounded(), "Not bounded");
        VX_ASSERT(isConvex(), "Not convex");
    }

    template<typename T, size_t DIM>
    std::vector<vx::HalfSpace<T, DIM>> vx::QuickHull<T, DIM>::getHalfSpaces() {
        std::vector<HalfSpace<T, DIM>>
                res;
        for (const auto &facet : facets)
            res.push_back(facet->half_space);
        return res;
    }

    template<typename T, size_t DIM>
    std::vector<std::array<size_t, DIM>> vx::QuickHull<T, DIM>::getFacets() {
        std::vector<std::array<size_t, DIM>> res;
        for (const auto &facet : facets)
            res.push_back(facet->vertices);
        return res;
    }

    template<typename T, size_t DIM>
    bool vx::QuickHull<T, DIM>::isConvex() const {
        // TODO: in general this does not guarantee an overall convexity. Only local convexity is checked
        for (const auto &f : facets) {
            for (size_t j = 0; j < DIM; j++) {
                if (f->neighbours[j]->isOutside(points[f->vertices[j]])) {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T, size_t DIM>
    bool vx::QuickHull<T, DIM>::isBounded() const {
        for (const auto &f : facets) {
            for (auto v: points) {
                if (f->isOutside(v)) {
                    std::cout << f->dist(v) << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

} // namespace vx

#endif //VORONOIX_ALG_CONVEX_HULL_KD_HPP
