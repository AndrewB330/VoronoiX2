#ifndef VORONOIX_ALG_DELAUNAY_2D_HPP
#define VORONOIX_ALG_DELAUNAY_2D_HPP

#include <vector>
#include <array>
#include <vx/geometry/vector.hpp>
#include <vx/geometry/polyhedron.hpp>
#include <vx/utils.hpp>
#include <algorithm>

#include "graph_delaunay.hpp"

namespace vx {

    template<typename T>
    class Delaunay2D {
    public:
        explicit Delaunay2D(const std::vector<vx::Vec<T, 2>> &points);

        const vx::DelaunayGraph_<T, 2> &getGraph() const;

    protected:
        void initialization();

        void finalization();

        size_t makeTriangle(const std::array<size_t, 3> &vertices);

        // fill in information about neighbour_facets to make this triangles connected
        inline void linkTriangles(size_t a_id, size_t b_id, size_t side_a, size_t side_b);

        // true - delaunay condition satisfied
        inline bool checkDelaunayCondition(size_t a_id, size_t b_id, size_t side_a, size_t side_b) const;

        // flip two neighbour triangles
        inline void flipTriangles(size_t a_id, size_t b_id, size_t side_a, size_t side_b);

        std::pair<size_t, size_t> queuePop();

        void queuePush(std::pair<size_t, size_t> v);

        // start BFS wave of flips
        void startFlipping();

        inline bool isBoundingPoint(size_t vertex) const;

        // there is 3 stages (default, aligning, finalizing)
        uint8_t flipping_stage = 0;

        // result facets of Delaunay triangulation
        std::vector<std::array<size_t, 3>> facets_vertices;
        std::vector<std::array<size_t, 3>> facets_neighbour;
        std::vector<std::array<size_t, 3>> facets_neighbour_sides;

        vx::DelaunayGraph_<T, 2> graph;

        // additional fields used to perform BFS on neighbours triangles
        std::vector<std::pair<size_t, size_t>> q;
        size_t q_tail = 0;
        size_t q_head = 0;

        // information about points used in computations
        size_t points_num = 0;
        std::vector<vx::Vec<T, 2>> points;

        // bounding box
        T min_x, max_x, min_y, max_y;
    };

    template<typename T>
    vx::Delaunay2D<T>::Delaunay2D(const std::vector<vx::Vec<T, 2>> &points_):points(points_) {
        VX_ASSERT(points.size() >= 3, "Number of points should be at least 3");

        initialization();

        size_t stripes_num = size_t(sqrt(0.2 * points.size()) + 1);
        T stripe_width = (max_x - min_x) / stripes_num + vx::epsilon<T>();
        T inverse_stripe_width = T(1) / stripe_width; // to speed-up computations

        std::vector<std::vector<size_t>> stripes(stripes_num);

        for (size_t i = 0; i < points_num; i++) {
            stripes[(size_t) ((points[i].x - min_x) * inverse_stripe_width)].push_back(i);
        }

        std::vector<size_t> order;
        for (size_t i = 0; i < stripes_num; i++) {
            std::sort(stripes[i].begin(), stripes[i].end(), [&](size_t a, size_t b) {
                return (i & 1u ? points[a].y > points[b].y : points[a].y < points[b].y);
            });
            order.insert(order.end(), stripes[i].begin(), stripes[i].end());
            std::vector<size_t>().swap(stripes[i]);
        }

        for (size_t vertex : order) {
            size_t triangle = facets_vertices.size() - 1; // last added triangle
            bool found = false;
            while (!found) {
                const auto &v = facets_vertices[triangle];
                if (cross_product(points[v[2]], points[vertex], points[v[1]]) < 0) {
                    triangle = facets_neighbour[triangle][0];
                } else if (cross_product(points[v[0]], points[vertex], points[v[2]]) < 0) {
                    triangle = facets_neighbour[triangle][1];
                } else if (cross_product(points[v[1]], points[vertex], points[v[0]]) < 0) {
                    triangle = facets_neighbour[triangle][2];
                } else {
                    found = true;
                }
                VX_ASSERT(triangle != vx::NO_PTR, "Vertex is outside of bounding box");
            }
            // inside
            size_t t0_id = triangle;
            size_t t1_id = makeTriangle({vertex, facets_vertices[triangle][2], facets_vertices[triangle][0]});
            size_t t2_id = makeTriangle({vertex, facets_vertices[triangle][0], facets_vertices[triangle][1]});

            facets_vertices[t0_id][0] = vertex;

            const auto &old_neighbours = facets_neighbour[t0_id];
            const auto &old_sides = facets_neighbour_sides[t0_id];

            linkTriangles(t1_id, old_neighbours[1], 0, old_sides[1]);
            linkTriangles(t2_id, old_neighbours[2], 0, old_sides[2]);

            linkTriangles(t0_id, t1_id, 1, 2);
            linkTriangles(t0_id, t2_id, 2, 1);
            linkTriangles(t1_id, t2_id, 1, 2);

            queuePush({t0_id, 0});
            queuePush({t1_id, 0});
            queuePush({t2_id, 0});

            startFlipping();
        }

        finalization();
    }

    template<typename T>
    bool Delaunay2D<T>::isBoundingPoint(size_t vertex) const {
        return vertex >= points_num;
    }

    template<typename T>
    size_t vx::Delaunay2D<T>::makeTriangle(const std::array<size_t, 3> &vertices) {
        facets_vertices.push_back(vertices);
        facets_neighbour.push_back({vx::NO_PTR, vx::NO_PTR, vx::NO_PTR});
        facets_neighbour_sides.push_back({0, 0, 0});
        return facets_vertices.size() - 1;
    }

    template<typename T>
    bool vx::Delaunay2D<T>::checkDelaunayCondition(size_t a_id, size_t b_id, size_t side_a, size_t side_b) const {
        if (a_id == vx::NO_PTR || b_id == vx::NO_PTR)
            return true;

        const auto &a = facets_vertices[a_id];
        const auto &b = facets_vertices[b_id];

        const auto &p0 = points[a[side_a]];
        const auto &p1 = points[a[vx::NEXT_PTR_3[side_a]]];
        const auto &p2 = points[a[vx::PREV_PTR_3[side_a]]];
        const auto &p3 = points[b[side_b]];

        if (flipping_stage == 1) {
            return !(isBoundingPoint(a[side_a]) && isBoundingPoint(a[vx::NEXT_PTR_3[side_a]]) &&
                     cross_product(p3, p2, p0) > 0);
        }
        if (flipping_stage == 2) {
            return !(!isBoundingPoint(a[side_a]) &&
                     isBoundingPoint(a[vx::NEXT_PTR_3[side_a]]) &&
                     !isBoundingPoint(a[vx::PREV_PTR_3[side_a]]) &&
                     cross_product(p3, p2, p0) > 0);
        }

        T d01x = p1.x - p0.x;
        T d01y = p1.y - p0.y;
        T d02x = p2.x - p0.x;
        T d02y = p2.y - p0.y;
        T d31x = p1.x - p3.x;
        T d31y = p1.y - p3.y;
        T d32x = p2.x - p3.x;
        T d32y = p2.y - p3.y;
        T cos_a = d32x * d31x + d32y * d31y;
        T cos_b = d01x * d02x + d01y * d02y;
        if (cos_a >= 0 && cos_b >= 0) // early check that both angles are less than 90
            return true;
        T sin_a = d32x * d31y - d32y * d31x;
        T sin_b = d01x * d02y - d01y * d02x;
        return sin_a * cos_b + sin_b * cos_a >= 0; // sum of angles less than 180
    }

    template<typename T>
    void Delaunay2D<T>::flipTriangles(size_t a_id, size_t b_id, size_t side_a, size_t side_b) {
        auto &a_vertices = facets_vertices[a_id];
        auto &b_vertices = facets_vertices[b_id];
        auto &a_neighbour = facets_neighbour[a_id];
        auto &b_neighbour = facets_neighbour[b_id];
        auto &a_neighbour_sides = facets_neighbour_sides[a_id];
        auto &b_neighbour_sides = facets_neighbour_sides[b_id];

        b_vertices[vx::PREV_PTR_3[side_b]] = a_vertices[side_a];
        a_vertices[vx::PREV_PTR_3[side_a]] = b_vertices[side_b];

        linkTriangles(a_id, b_neighbour[vx::NEXT_PTR_3[side_b]], side_a,
                      b_neighbour_sides[vx::NEXT_PTR_3[side_b]]);
        linkTriangles(b_id, a_neighbour[vx::NEXT_PTR_3[side_a]], side_b,
                      a_neighbour_sides[vx::NEXT_PTR_3[side_a]]);

        linkTriangles(a_id, b_id, vx::NEXT_PTR_3[side_a], vx::NEXT_PTR_3[side_b]);
    }

    template<typename T>
    void Delaunay2D<T>::initialization() {
        points_num = points.size();
        size_t max_number_of_triangles = 2 * points_num + 8;

        facets_vertices.reserve(max_number_of_triangles);
        facets_neighbour.reserve(max_number_of_triangles);
        facets_neighbour_sides.reserve(max_number_of_triangles);

        q.reserve(3*max_number_of_triangles + 100); // todo: not enough in WORST case

        min_x = points[0].x, max_x = points[0].x, min_y = points[0].y, max_y = points[0].y;
        for (const auto &p : points) {
            min_x = std::min(min_x, p.x - vx::epsilon<T>());
            max_x = std::max(max_x, p.x + vx::epsilon<T>());
            min_y = std::min(min_y, p.y - vx::epsilon<T>());
            max_y = std::max(max_y, p.y + vx::epsilon<T>());
        }

        T width = (max_x - min_x);
        T height = (max_y - min_y);

        points.emplace_back(min_x - 0.07 * width, min_y - 0.07 * height); // left bottom
        points.emplace_back(max_x + 0.06 * width, min_y - 0.06 * height); // right bottom
        points.emplace_back(max_x + 0.08 * width, max_y + 0.08 * height); // right top
        points.emplace_back(min_x - 0.06 * width, max_y + 0.05 * height); // left top

        makeTriangle({points_num + 0, points_num + 1, points_num + 2}); // bottom bottom corner
        makeTriangle({points_num + 2, points_num + 3, points_num + 0}); // left top corner
        linkTriangles(0, 1, 1, 1);
        if (!checkDelaunayCondition(0, 1, 1, 1))
            flipTriangles(0, 1, 1, 1);
    }

    template<typename T>
    void Delaunay2D<T>::linkTriangles(size_t a_id, size_t b_id, size_t side_a, size_t side_b) {
        if (a_id != vx::NO_PTR) {
            facets_neighbour[a_id][side_a] = b_id;
            facets_neighbour_sides[a_id][side_a] = side_b;
        }
        if (b_id != vx::NO_PTR) {
            facets_neighbour[b_id][side_b] = a_id;
            facets_neighbour_sides[b_id][side_b] = side_a;
        }
    }

    template<typename T>
    void Delaunay2D<T>::startFlipping() {
        while (q_head < q_tail) {
            auto[cur, side] = queuePop();
            size_t neighbour = facets_neighbour[cur][side];
            size_t neighbour_side = facets_neighbour_sides[cur][side];
            if (neighbour != vx::NO_PTR) {
                if (!checkDelaunayCondition(cur, neighbour, side, neighbour_side)) {
                    flipTriangles(cur, neighbour, side, neighbour_side);
                    queuePush({neighbour, vx::PREV_PTR_3[neighbour_side]});
                    queuePush({cur, side});
                }
            }
        }
        q_head = q_tail = 0;
    }

    template<typename T>
    void Delaunay2D<T>::finalization() {
        flipping_stage = 1;
        for (size_t i = 0; i < facets_vertices.size(); i++) {
            for (size_t side = 0; side < 3; side++) {
                if (isBoundingPoint(facets_vertices[i][side]) &&
                    isBoundingPoint(facets_vertices[i][vx::NEXT_PTR_3[side]])) {
                    queuePush({i, side});
                }
            }
        }
        startFlipping();

        flipping_stage = 2;
        for (size_t i = 0; i < facets_vertices.size(); i++) {
            for (size_t side = 0; side < 3; side++) {
                if (!isBoundingPoint(facets_vertices[i][side]) &&
                    isBoundingPoint(facets_vertices[i][vx::NEXT_PTR_3[side]]) &&
                    !isBoundingPoint(facets_vertices[i][vx::PREV_PTR_3[side]])) {
                    queuePush({i, side});
                }
            }
        }
        startFlipping();

        std::vector<size_t> new_id(facets_vertices.size());
        size_t valid_faces_num = 0;
        std::vector<size_t> bb;
        for (size_t i = 0; i < facets_vertices.size(); i++) {
            bool is_valid = std::all_of(
                    facets_vertices[i].begin(),
                    facets_vertices[i].end(),
                    [&](auto v) { return !isBoundingPoint(v); }
            );
            if (is_valid) {
                facets_vertices[valid_faces_num] = facets_vertices[i];
                facets_neighbour[valid_faces_num] = facets_neighbour[i];
                facets_neighbour_sides[valid_faces_num] = facets_neighbour_sides[i];
                new_id[i] = valid_faces_num;
                ++valid_faces_num;
            } else {
                new_id[i] = vx::NO_PTR;
                for(int k = 0; k < 3; k++)
                if (!isBoundingPoint(facets_vertices[i][k]))
                    bb.push_back(facets_vertices[i][k]);
            }
        }
        std::sort(bb.begin(), bb.end());
        auto en = std::unique(bb.begin(), bb.end());
        bb.erase(en, bb.end());
        std::vector<Vec<T, 2>> pts;
        for(auto i : bb) pts.push_back(points[i]);
        vx::GrahamScan<T> ch(pts);
        std::cerr << ch.getVertices().size() << std::endl;

        facets_vertices.resize(valid_faces_num);
        facets_neighbour.resize(valid_faces_num);
        facets_neighbour_sides.resize(valid_faces_num);
        for (auto &f : facets_neighbour) {
            for (auto &neighbour : f) {
                neighbour = new_id[neighbour];
            }
        }

        graph.facets_vertices.swap(facets_vertices);
        graph.facets_adjacent.swap(facets_neighbour);
        graph.facets_adjacent_side.swap(facets_neighbour_sides);
        graph.points.swap(points);

        graph.adjacent_points.resize(graph.points.size());
        for(int i = 0; i < graph.points.size(); i++) {
            graph.adjacent_points[i].reserve(4);
        }

        for(const auto & vertices : graph.facets_vertices) {
            for(uint8_t side = 0; side < 3; side++) {
                graph.adjacent_points[vertices[side]].push_back(vertices[vx::NEXT_PTR_3[side]]);
                graph.adjacent_points[vertices[side]].push_back(vertices[vx::PREV_PTR_3[side]]);
            }
        }
        // todo: vertices and edges
    }

    template<typename T>
    std::pair<size_t, size_t> Delaunay2D<T>::queuePop() {
        return q[q_head++];
    }

    template<typename T>
    void Delaunay2D<T>::queuePush(std::pair<size_t, size_t> v) {
        q[q_tail++] = v;
    }

    template<typename T>
    const vx::DelaunayGraph_<T, 2> &Delaunay2D<T>::getGraph() const {
        return graph;
    }

    template<typename T>
    std::vector<vx::Polyhedron<T, 2>> voronoi(const std::vector<Vec<T, 2>> & points, T offset) {
        std::vector<vx::Polyhedron<T, 2>> res;
        auto delaunay = Delaunay2D<T>(points).getGraph();
        for(int i = 0; i < delaunay.numVertices(); i++) {
            std::vector<HalfSpace<T, 2>> hs;
            hs.emplace_back(vx::Vec<T, 2>(1, 0), 4000);
            hs.emplace_back(vx::Vec<T, 2>(-1, 0), 4000);
            hs.emplace_back(vx::Vec<T, 2>(0, 1), 4000);
            hs.emplace_back(vx::Vec<T, 2>(0, -1), 4000);
            auto a = delaunay.getVertex(i).getPoint();
            for(DelaunayVertex<T, 2> & v : delaunay.getVertex(i).getAdjacentVertices()) {
                auto b = v.getPoint();
                auto dif = b - a;
                auto len = length(dif);
                auto dir = dif / len;
                hs.emplace_back(dir, len/2 - offset + dot(a, dir));
            }
            res.push_back(makePolyhedron(hs, a));
        }
        return res;
    }

} // namespace vx

#endif //VORONOIX_ALG_DELAUNAY_2D_HPP
