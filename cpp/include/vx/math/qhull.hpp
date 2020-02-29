#ifndef VORONOIX_QHULL_HPP
#define VORONOIX_QHULL_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <ctime>
#include "vec.hpp"
#include "geometry.hpp"

namespace vx {

    const int ulp = 1024 * 128;

    template<typename T, size_t DIM, size_t CNT>
    struct BiggestSimplexBuilder {
        const std::vector<Vec<T, DIM>> &pts;

        explicit BiggestSimplexBuilder(const std::vector<Vec<T, DIM>> &pts) : pts(pts) {}

        std::set<size_t> get_biggest_simplex() {
            BiggestSimplexBuilder<T, DIM, CNT - 1> simplex_builder(pts);
            auto res = simplex_builder.get_biggest_simplex();
            std::vector<size_t> res_vec(res.begin(), res.end());
            T max_volume = T(0);
            size_t max_volume_index = pts.size();
            for (size_t i = 0; i < pts.size(); i++) {
                if (res.find(i) != res.end()) continue;
                vx::Mat<T, CNT - 1, CNT - 1> gramian_matrix;
                for (size_t ii = 0; ii < CNT - 1; ii++) {
                    for (size_t jj = 0; jj < CNT - 1; jj++) {
                        gramian_matrix[ii][jj] = dot(pts[res_vec[ii]] - pts[i], pts[res_vec[jj]] - pts[i]);
                    }
                }
                double volume = det(gramian_matrix);
                if (volume > max_volume) {
                    max_volume = volume;
                    max_volume_index = i;
                }
            }
            res.insert(max_volume_index);
            return res;
        }
    };

    template<typename T, size_t DIM>
    struct BiggestSimplexBuilder<T, DIM, 1> {
        const std::vector<Vec<T, DIM>> &pts;

        explicit BiggestSimplexBuilder(const std::vector<Vec<T, DIM>> &pts) : pts(pts) {}

        std::set<size_t> get_biggest_simplex() {
            T max_dist = T(0);
            size_t max_dist_index = 0;
            Vec<T, DIM> c;
            for (auto p: pts) c += p;
            c /= pts.size();
            for (size_t i = 0; i < pts.size(); i++) {
                T dist = dist_square(c, pts[i]);
                if (dist > max_dist) {
                    max_dist = dist;
                    max_dist_index = i;
                }
            }
            std::set<size_t> res;
            res.insert(max_dist_index);
            return res;
        }
    };

    template<typename T, size_t DIM>
    struct Facet {
        using FacetPtr = std::shared_ptr<Facet<T, DIM>>;
        size_t vertices[DIM] = {};
        FacetPtr neighbours[DIM] = {};

        size_t furthest = 0;
        T furthest_dist = T(0);
        std::vector<size_t> outside;
        std::vector<size_t> coplanar;

        HalfSpace<T, DIM> half_space;

        explicit Facet() {}

        void replaceNeighbour(FacetPtr old, FacetPtr current) {
            auto it = std::find(neighbours, neighbours + DIM, old);
            VX_ASSERT(it != neighbours + DIM, "Whops");
            *it = current;
        }

        void updateHalfSpace(const std::vector<Vec<T, DIM>> &p, const GeometryContext<T> &ctx) {
            std::vector<Vec<T, DIM>> basis(DIM);
            for (size_t i = 0; i < DIM; i++)
                basis[i] = p.at(vertices[i]);
            half_space = makeHalfSpace(basis, ctx);
        }

        void checkAndPutVertex(size_t i, const Vec<T, DIM> &v) {
            if (isCoplanar(v)) {
                //coplanar.push_back(i);
            } else if (isOutside(v)) {
                outside.push_back(i);
                T cur_dist = dist(v);
                if (cur_dist > furthest_dist) {
                    furthest_dist = cur_dist;
                    furthest = i;
                }
            }
        }

        void flip() {
            std::swap(vertices[0], vertices[1]);
            std::swap(neighbours[0], neighbours[1]);
            half_space.flip();
        }

        T dist(const Vec<T, DIM> &v) const {
            return half_space.orientedDistance(v);
        }

        bool isOutside(const Vec<T, DIM> &v) const {
            return half_space.isAbove(v);
        }

        bool isCoplanar(const Vec<T, DIM> &v) const {
            return half_space.isOn(v);
        }
    };

    template<typename T, size_t DIM>
    struct Ridge {
        size_t vertices[DIM - 1];

        Ridge(const std::shared_ptr<Facet<T, DIM>> facet, size_t side) {
            auto ptr = vertices;
            for (size_t i = 0; i < DIM; i++) {
                if (i != side) {
                    *(ptr++) = facet->vertices[i];
                }
            }
            std::sort(vertices, vertices + DIM - 1);
        }

        bool operator==(const Ridge &r) const {
            return std::equal(vertices, vertices + DIM - 1, r.vertices);
        }
    };


    template<typename T, size_t DIM>
    struct RidgeHasher {
        size_t operator()(const Ridge<T, DIM> &r) const {
            size_t h = 0;
            for (size_t i = 0; i < DIM - 1; i++) {
                h ^= (r.vertices[i] * 331) ^ (r.vertices[i]);
            }
            return h;
        }
    };

    template<typename T, size_t DIM>
    bool operator<(const Facet<T, DIM> &a, const Facet<T, DIM> &b) {
        return a.furthest_dist > b.furthest_dist;
    }

    template<typename T, size_t DIM>
    struct QHull {

        Vec<T, DIM> inside;
        GeometryContext<T> ctx;
        std::vector<Vec<T, DIM>> points;
        std::set<std::shared_ptr<Facet<T, DIM>>> facet_candidates;
        std::unordered_set<std::shared_ptr<Facet<T, DIM>>> facets;


        void addFacet(const std::shared_ptr<Facet<T, DIM>> &facet) {
            facets.insert(facet);
            if (!facet->outside.empty())
                facet_candidates.insert(facet);
        }

        void deleteFacet(const std::shared_ptr<Facet<T, DIM>> &facet) {
            facet_candidates.erase(facet);
            facets.erase(facet);
            for (size_t i = 0; i < DIM; i++) {
                facet->neighbours[i] = nullptr;
            }
        }

        explicit QHull(const std::vector<vx::Vec<T, DIM>> &points) : points(points) {
            std::cout << std::clock() * 1.0 / CLOCKS_PER_SEC << std::endl;
            VX_ASSERT(points.size() > DIM, "Number of points should be at least DIM+1 for QHull algrithm");
            buildInitialFacets();

            std::unordered_set<size_t> outside_vertex_candidates{};
            while (!facet_candidates.empty()) {
                outside_vertex_candidates.clear();
                auto facet = *facet_candidates.begin();
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

                    std::copy_n(old_facet->vertices, DIM, new_facet->vertices);
                    new_facet->vertices[connection] = furthest_vertex;
                    new_facet->neighbours[connection] = old_facet->neighbours[connection];
                    new_facet->neighbours[connection]->replaceNeighbour(old_facet, new_facet);

                    new_facet->updateHalfSpace(points, ctx);

                    new_facet->outside.reserve(outside_vertex_candidates.size() / 2);
                    for (size_t i : outside_vertex_candidates)
                        new_facet->checkAndPutVertex(i, points[i]);
                    // TODO: should we check only candidates?not all?

                    new_facets.push_back(new_facet);
                    addFacet(new_facet);
                    VX_ASSERT(!new_facet->isOutside(inside), "New facet should be directed outside");
                }

                connectAllNewFacets(new_facets);

                for (auto to_delete: visible_set) deleteFacet(to_delete);

                //std::cout << facets.size() << std::endl;
                /*for (auto f: facets) {
                    for (size_t j = 0; j < DIM; j++) {
                        auto n = f->neighbours[j];
                        VX_ASSERT(n != nullptr, "No neighbour");
                        auto it = std::find(n->neighbours, n->neighbours + DIM, f);
                        VX_ASSERT(it != n->neighbours + DIM, "Neighbours not consistent");
                    }
                }*/
                //VX_ASSERT(isConvex(), "Not convex");
            }
            std::cout << std::clock() * 1.0 / CLOCKS_PER_SEC << std::endl;
            VX_ASSERT(isBounded(), "Not bounded");
            VX_ASSERT(isConvex(), "Not convex");
        }

        void connectAllNewFacets(const std::vector<std::shared_ptr<Facet<T, DIM>>> &new_facets) {
            using RidgeInfo = std::pair<std::shared_ptr<Facet<T, DIM>>, size_t>;
            std::unordered_map<Ridge<T, DIM>, RidgeInfo, RidgeHasher<T, DIM>> ridge_map{};
            for (auto facet : new_facets) {
                for (size_t side = 0; side < DIM; side++) {
                    auto ridge = Ridge<T, DIM>(facet, side);
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

        void buildInitialFacets() {
            BiggestSimplexBuilder<T, DIM, DIM + 1> biggest_simplex(points);
            auto simplex = biggest_simplex.get_biggest_simplex();
            std::vector<size_t> simplex_vec(simplex.begin(), simplex.end());

            for (auto v : simplex)
                inside += points[v];
            inside /= simplex.size();

            std::vector<std::shared_ptr<Facet<T, DIM>>> pool(simplex_vec.size());
            for (size_t i = 0; i < simplex_vec.size(); i++) {
                pool[i] = std::make_shared<Facet<T, DIM>>();
            }
            for (size_t k = 0; k < simplex_vec.size(); k++) {
                auto facet = pool[k];
                for (size_t i = 0; i < DIM + 1; i++) {
                    if (i != k) {
                        size_t ii = i - (i > k);
                        facet->vertices[ii] = simplex_vec[i];
                        facet->neighbours[ii] = pool[i];
                    }
                }
                facet->updateHalfSpace(points, ctx);
                if (facet->isOutside(inside)) {
                    facet->flip();
                }
                for (size_t i = 0; i < points.size(); i++)
                    facet->checkAndPutVertex(i, points[i]);
                addFacet(facet);
            }
        }

        bool isConvex() const {
            // TODO: in general this does not guarantee an overall convexity. Only local convexity is checked
            for (auto f : facets) {
                for (size_t j = 0; j < DIM; j++) {
                    if (f->neighbours[j]->isOutside(points[f->vertices[j]])) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool isBounded() const {
            for (auto f : facets) {
                for (auto v: points) {
                    if (f->isOutside(v)) {
                        return false;
                    }
                }
            }
            return true;
        }

        void save_json(const std::string &filename) {
            std::ofstream out(filename);
            if (DIM == 4) {
                out << "let data = { \"v\":[\n";
                for (auto p : points) {
                    out << "[" << p[0] << "," << p[1] << "," << p[2] << "]," << std::endl;
                }
                out << "], \"f\": [\n";
                for (auto f : facets) {
                    if (f->half_space.normal[DIM - 1] < 0) {
                        out << "[" << f->vertices[0] << ", " << f->vertices[1] << ", " << f->vertices[2] << "],"
                            << std::endl;
                        out << "[" << f->vertices[0] << ", " << f->vertices[1] << ", " << f->vertices[3] << "],"
                            << std::endl;
                        out << "[" << f->vertices[0] << ", " << f->vertices[2] << ", " << f->vertices[3] << "],"
                            << std::endl;
                        out << "[" << f->vertices[1] << ", " << f->vertices[2] << ", " << f->vertices[3] << "],"
                            << std::endl;
                    }
                }
                out << "]}";
            }
            if (DIM == 3) {
                out << "let data = { \"v\":[\n";
                for (auto p : points) {
                    out << "[" << p[0] << "," << p[1] << "," << p[2] << "]," << std::endl;
                }
                out << "[0,0,0]], \"f\": [\n";
                for (auto f : facets) {
                    //if (f.half_space.normal[DIM - 1] < 0)
                    out << "[" << f->vertices[0] << ", " << f->vertices[1] << ", " << f->vertices[2]
                        << "],"
                        << std::endl;
                }
                out << "]}";
            }

        }

        Polyhedron<T, DIM> getHull() {
            Polyhedron<T, DIM> res;
            res.points = points;
            res.ctx = ctx;
            for (auto facet: facets) {
                VerticesList<DIM> vertices{};
                std::copy_n(facet->vertices, DIM, vertices.vertices);
                res.facets.push_back(vertices);
            }
            return res;
        }

    };


    template<typename T, size_t DIM>
    Polyhedron<T, DIM> makePolyhedron(const std::vector<HalfSpace<T, DIM>> &half_spaces, const Vec<T, DIM> &center) {
        std::vector<Vec<T, DIM>> dual_vertices;
        for (auto hs : half_spaces) {
            auto offset = hs.offset - dot(center, hs.normal);
            dual_vertices.push_back(hs.normal / offset);
        }
        QHull<T, DIM> dual_qhull(dual_vertices);
        std::vector<Vec<T, DIM>> result_vertices;
        for (auto facet : dual_qhull.facets) {
            result_vertices.push_back(facet->half_space.normal / facet->half_space.offset + center);
        }
        QHull<T, DIM> result_qhull(result_vertices);
        return result_qhull.getHull();
    }


}


#endif //VORONOIX_QHULL_HPP
