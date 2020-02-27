#ifndef VORONOIX_QHULL_HPP
#define VORONOIX_QHULL_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include "vec.hpp"

namespace vx {

    const int ulp = 1024 * 128;

    template<typename T, size_t DIM>
    struct Context {
        size_t n; // number of points
        const std::vector<vx::Vec<T, DIM>> &p; // points
        vx::Vec<T, DIM> c; // centroid
        T eps = std::numeric_limits<T>::epsilon() * ulp;
        T tigh_eps = std::numeric_limits<T>::epsilon() * (ulp / 4);
        T span = T(1);

        explicit Context(const std::vector<vx::Vec<T, DIM>> &points) : p(points) {
            for (const auto &v : points) c += v;
            c /= points.size();
            n = points.size();
            span = T(0);
            for (const auto &v : points) span = std::max(span, 2 * dist(c, v));
            eps = span * (std::numeric_limits<T>::epsilon() * ulp);
            tigh_eps = span * (std::numeric_limits<T>::epsilon() * (ulp / 4));
        }
    };

    template<typename T, size_t DIM, size_t CNT>
    struct BiggestSimplexBuilder {
        const Context<T, DIM> &ctx;

        explicit BiggestSimplexBuilder(const Context<T, DIM> &ctx) : ctx(ctx) {}

        std::set<size_t> get_biggest_simplex() {
            BiggestSimplexBuilder<T, DIM, CNT - 1> simplex_builder(ctx);
            auto res = simplex_builder.get_biggest_simplex();
            std::vector<size_t> res_vec(res.begin(), res.end());
            T max_volume = T(0);
            size_t max_volume_index = ctx.p.size();
            for (size_t i = 0; i < ctx.n; i++) {
                if (res.find(i) != res.end()) continue;
                vx::Mat<T, CNT - 1, CNT - 1> gramian_matrix;
                for (size_t ii = 0; ii < CNT - 1; ii++) {
                    for (size_t jj = 0; jj < CNT - 1; jj++) {
                        gramian_matrix[ii][jj] = dot(ctx.p[res_vec[ii]] - ctx.p[i], ctx.p[res_vec[jj]] - ctx.p[i]);
                    }
                }
                double volume = det(gramian_matrix);
                /*if (DIM + 1 == CNT) {
                    Mat<T, DIM, DIM> basis;
                    for(size_t i = 0; i < DIM; i++) basis.set_row(i, ctx.p[res_vec[i]] - ctx.p[i]);
                    volume = sdet(basis);
                }*/
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
        const Context<T, DIM> &ctx;

        explicit BiggestSimplexBuilder(const Context<T, DIM> &ctx) : ctx(ctx) {}

        std::set<size_t> get_biggest_simplex() {
            T max_dist = T(0);
            size_t max_dist_index = 0;
            for (size_t i = 0; i < ctx.n; i++) {
                T dist = dist_square(ctx.c, ctx.p[i]);
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
        const Context<T, DIM> *ctx;
        size_t vertices[DIM] = {};
        size_t neighbours[DIM] = {};

        size_t furthest = 0;
        T furthest_dist = T(0);
        std::vector<size_t> outside;
        std::vector<size_t> coplanar;

        Vec <T, DIM> normal;
        T offset = T(0);

        explicit Facet(const Context<T, DIM> &ctx) : ctx(&ctx) {}

        void update_equation() {
            Mat<T, DIM, DIM> mat;
            for (size_t i = 0; i < DIM; i++) {
                mat.set_col(i, ctx->p[vertices[i]]);
            }
            offset = det(mat);
            Mat<T, DIM, DIM> mat_copy;
            for (size_t i = 0; i < DIM; i++) {
                mat_copy = mat;
                mat_copy.set_row(i, vec_ones<T, DIM>());
                normal[i] = det(mat_copy);
            }
            T d = length(normal);
            normal /= d;
            offset /= d;
            for (size_t v : vertices) {

                if (!(dot(ctx->p[v], normal) < offset + ctx->tigh_eps &&
                      dot(ctx->p[v], normal) > offset - ctx->tigh_eps &&
                      "Face vertices must satisfy plane equation")) {
                    /*std::cout << ctx->p[vertices[0]] << std::endl;
                    std::cout << ctx->p[vertices[1]] << std::endl;
                    std::cout << "ohh" << std::endl;*/
                }
            }
        }

        void update_outer_vertices() {
            for (size_t i = 0; i < ctx->n; i++) {
                if (check_coplanar(i)) {
                    coplanar.push_back(i);
                } else if (check_outside(i)) {
                    outside.push_back(i);
                    T cur_dist = dist(i);
                    if (cur_dist > furthest_dist) {
                        furthest_dist = cur_dist;
                        furthest = i;
                    }
                }
            }
        }

        void update_outer_vertices(const Facet &parent) {
            for (size_t i = 0; i < ctx->n; i++) {
                if (check_coplanar(i)) {
                    coplanar.push_back(i);
                } else if (check_outside(i)) {
                    outside.push_back(i);
                    T cur_dist = dist(i);
                    if (cur_dist > furthest_dist) {
                        furthest_dist = cur_dist;
                        furthest = i;
                    }
                }
            }
        }

        void flip() {
            std::swap(vertices[0], vertices[1]);
            std::swap(neighbours[0], neighbours[1]);
            normal = -normal;
            offset = -offset;
        }

        T dist(size_t vertex) const {
            return vx::dot(ctx->p[vertex], normal) - offset;
        }

        bool check_outside(size_t vertex) const {
            return dot(ctx->p[vertex], normal) > offset + ctx->eps;
        }

        bool check_outside(const Vec <T, DIM> &v) const {
            return dot(v, normal) > offset + ctx->eps;
        }

        bool check_coplanar(size_t vertex) const {
            return dot(ctx->p[vertex], normal) <= offset + ctx->eps &&
                   dot(ctx->p[vertex], normal) >= offset - ctx->eps;
        }

        bool check_coplanar(const Vec <T, DIM> &v) const {
            return dot(v, normal) <= offset + ctx->eps &&
                   dot(v, normal) >= offset - ctx->eps;
        }
    };

    template<typename T, size_t DIM>
    struct QHull {

        Context<T, DIM> ctx;
        std::vector<Facet<T, DIM>> facets;

        explicit QHull(const std::vector<vx::Vec<T, DIM>> &points) : ctx(points) {
            facets = get_initial_facets();
            std::cout << ctx.p.size() << std::endl;
            for (size_t facet = choose_facet(); facet < facets.size(); facet = choose_facet()) {
                size_t furthest_vertex = facets[facet].furthest;
                std::vector<size_t> visible_set;
                std::vector<char> visited(facets.size(), 0);
                std::queue<size_t> queue;
                queue.push(facet);
                visited[facet] = true;
                std::vector<Facet<T, DIM>> new_facets;
                std::vector<size_t *> new_facets_neighbours;

                while (!queue.empty()) {
                    size_t current_facet = queue.front();
                    queue.pop();
                    visible_set.push_back(current_facet);
                    for (size_t i = 0; i < DIM; i++) {
                        size_t neighbour = facets[current_facet].neighbours[i];
                        if (!facets[neighbour].check_outside(furthest_vertex)) {
                            Facet ridge = facets[current_facet];
                            //std::fill(ridge.neighbours, ridge.neighbours + DIM, facets.size());
                            ridge.vertices[i] = furthest_vertex;
                            ridge.neighbours[i] = neighbour;
                            ridge.outside.clear();
                            ridge.coplanar.clear();
                            ridge.furthest_dist = T(0);
                            new_facets.push_back(ridge);
                            size_t index = std::find(facets[neighbour].neighbours,
                                                     facets[neighbour].neighbours + DIM,
                                                     current_facet)
                                           - facets[neighbour].neighbours;
                            new_facets_neighbours.push_back(&facets[neighbour].neighbours[index]);
                            if (ridge.check_outside(ctx.c)) {
                                std::cout << ctx.p.size() << std::endl;
                            }
                        } else if (!visited[neighbour]) {
                            queue.push(neighbour);
                            visited[neighbour] = true;
                        }
                    }
                }

                std::vector<size_t> new_ids = visible_set;
                for (size_t i = 0; i + visible_set.size() < new_facets.size(); i++)
                    new_ids.push_back(facets.size() + i);

                for (size_t cur_facet = 0; cur_facet < new_facets.size(); cur_facet++) {
                    *new_facets_neighbours[cur_facet] = new_ids[cur_facet];
                    for (size_t other_facet = 0; other_facet < new_facets.size(); other_facet++) {
                        if (other_facet == cur_facet) continue;
                        std::vector<size_t> tmp_vertices(new_facets[cur_facet].vertices,
                                                         new_facets[cur_facet].vertices + DIM);
                        std::vector<size_t> tmp_vertices_2(new_facets[other_facet].vertices,
                                                           new_facets[other_facet].vertices + DIM);
                        std::sort(tmp_vertices.begin(),
                                  tmp_vertices.end());
                        std::sort(tmp_vertices_2.begin(),
                                  tmp_vertices_2.end());
                        auto end = std::set_difference(tmp_vertices.begin(),
                                                       tmp_vertices.end(),
                                                       tmp_vertices_2.begin(),
                                                       tmp_vertices_2.end(),
                                                       tmp_vertices.begin());
                        tmp_vertices.erase(end, tmp_vertices.end());
                        if (tmp_vertices.size() == 1) {
                            size_t index = std::find(new_facets[cur_facet].vertices,
                                                     new_facets[cur_facet].vertices + DIM,
                                                     tmp_vertices[0]) - new_facets[cur_facet].vertices;
                            new_facets[cur_facet].neighbours[index] = new_ids[other_facet];
                        }
                    }
                    new_facets[cur_facet].update_equation();
                    for (size_t visible_facet : visible_set) {
                        new_facets[cur_facet].update_outer_vertices(facets[visible_facet]);
                    }
                }
                for (size_t i = 0; i < new_facets.size(); i++) {
                    if (new_ids[i] < facets.size())
                        facets[new_ids[i]] = new_facets[i];
                    else
                        facets.push_back(new_facets[i]);
                }
            }
            assert(check_hull());
            assert(check_convex());
        }

        size_t choose_facet() {
            T best_facet_dist = 0;
            size_t best_facet = facets.size();
            for (size_t i = 0; i < facets.size(); i++) {
                if (!facets[i].outside.empty() && facets[i].furthest_dist > best_facet_dist) {
                    best_facet_dist = facets[i].furthest_dist;
                    best_facet = i;
                }
            }
            return best_facet;
        }

        std::vector<Facet<T, DIM>> get_initial_facets() {
            std::vector<Facet<T, DIM>> res;
            BiggestSimplexBuilder<T, DIM, DIM + 1> biggest_simplex(ctx);
            auto simplex = biggest_simplex.get_biggest_simplex();
            std::vector<size_t> simplex_vec(simplex.begin(), simplex.end());

            Vec<T, DIM> centroid;
            for (auto v : simplex) centroid += ctx.p[v];
            centroid /= simplex.size();
            ctx.c = centroid;

            for (size_t k = 0; k < simplex_vec.size(); k++) {
                Facet<T, DIM> facet(ctx);
                for (size_t i = 0; i < DIM + 1; i++) {
                    if (i != k) {
                        size_t ii = i - (i > k);
                        facet.vertices[ii] = simplex_vec[i];
                        facet.neighbours[ii] = i;
                    }
                }
                facet.update_equation();
                if (facet.check_outside(ctx.c)) {
                    facet.flip();
                }
                facet.update_outer_vertices();
                res.push_back(facet);
            }
            return res;
        }

        bool check_convex() const {
            for (size_t i = 0; i < facets.size(); i++) {
                for (size_t j = 0; j < DIM; j++) {
                    if (facets[facets[i].neighbours[j]].check_outside(facets[i].vertices[j])) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool check_hull() const {
            for (size_t facet = 0; facet < facets.size(); facet++) {
                for (size_t i = 0; i < ctx.n; i++) {
                    if (facets[facet].check_outside(i))
                        return false;
                }
            }
            return true;
        }

        void save_json(const std::string &filename) {
            std::ofstream out(filename);
            if (DIM == 4) {
                out << "let data = { \"v\":[\n";
                for (auto p : ctx.p) {
                    out << "[" << p[0] << "," << p[1] << "," << p[2] << "]," << std::endl;
                }
                out << "], \"f\": [\n";
                for (auto f : facets) {
                    if (f.normal[DIM - 1] < 0)
                    {
                        out << "[" << f.vertices[0] << ", " << f.vertices[1] << ", " << f.vertices[2] << "],"
                            << std::endl;
                        out << "[" << f.vertices[0] << ", " << f.vertices[1] << ", " << f.vertices[3] << "],"
                            << std::endl;
                        out << "[" << f.vertices[0] << ", " << f.vertices[2] << ", " << f.vertices[3] << "],"
                            << std::endl;
                        out << "[" << f.vertices[1] << ", " << f.vertices[2] << ", " << f.vertices[3] << "],"
                            << std::endl;
                    }
                }
                out << "]}";
            }
            if (DIM == 3) {
                out << "let data = { \"v\":[\n";
                for (auto p : ctx.p) {
                    out << "[" << p[0] << "," << p[1] << "," << p[2] << "]," << std::endl;
                }
                out << "[0,0,0]], \"f\": [\n";
                for (size_t i = 0; i < facets.size(); i++) {
                    auto f = facets[i];
                    //if (f.normal[DIM - 1] < 0)
                        out << "[" << f.vertices[0] << ", " << f.vertices[1] << ", " << f.vertices[2]
                            << (i + 1 == facets.size() ? "]": "],")
                            << std::endl;
                }
                out << "]}";
            }

        }

    };


}


#endif //VORONOIX_QHULL_HPP
