#ifndef VORONOIX_VORONOID_2D_HPP
#define VORONOIX_VORONOID_2D_HPP

#include "alg_delaunay_2d.hpp"
#include "graph_voronoi.hpp"

namespace vx {

    template<typename T>
    class Voronoi2D {
    public:
        explicit Voronoi2D(const std::vector<vx::Vec<T, 2>> &points);

        const vx::VoronoiGraph<T, 2> &getGraph() const;

    public:
        vx::Vec<T, 2>
        getCircumcircleCenter(const std::array<size_t, 3> &vertices, const std::vector<vx::Vec<T, 2>> &points);

        bool buildRegion(size_t site, size_t triangle, size_t side, const vx::DelaunayGraph_<T, 2> &delaunay);

        vx::VoronoiGraph<T, 2> graph;
    };

    template<typename T>
    Voronoi2D<T>::Voronoi2D(const std::vector<vx::Vec<T, 2>> &points) {
        vx::Delaunay2D<T> delaunay(points);

        auto delaunay_graph = delaunay.getGraph();//todo: move

        auto &junction_points = graph.junctions;
        auto &adjacent_sites = graph.adjacent_sites;
        auto &adjacent_ridges = graph.adjacent_sites;
        auto &ridges = graph.ridges;

        junction_points.reserve(delaunay_graph.numFacets());
        for (size_t i = 0; i < delaunay_graph.numFacets(); i++) {
            junction_points.push_back(getCircumcircleCenter(delaunay_graph.facets_vertices[i], points));
        }

        adjacent_sites.resize(points.size());
        adjacent_ridges.resize(points.size());

        for (size_t i = 0; i < delaunay_graph.numFacets(); i++) {
            for (size_t side = 0; side < 3; side++) {
                buildRegion(delaunay_graph.facets_vertices[i][side], i, side, delaunay_graph);
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            if (adjacent_ridges[i].size() == 1)
                adjacent_ridges[i].clear();
        }

        graph.points.swap(delaunay_graph.points);
    }

    template<typename T>
    vx::Vec<T, 2>
    Voronoi2D<T>::getCircumcircleCenter(const std::array<size_t, 3> &vertices,
                                        const std::vector<vx::Vec<T, 2>> &points) {
        const auto &A = points[vertices[0]];
        const auto &B = points[vertices[1]];
        const auto &C = points[vertices[2]];
        T area = std::abs(cross_product(A, B, C)) / 2;
        T a = dist_square(B, C) / (8 * area * area) * dot(A - B, A - C);
        T b = dist_square(C, A) / (8 * area * area) * dot(B - A, B - C);
        T c = dist_square(A, B) / (8 * area * area) * dot(C - A, C - B);
        return A * a + B * b + C * c;
    }

    template<typename T>
    bool
    Voronoi2D<T>::buildRegion(size_t site, size_t triangle, size_t side, const vx::DelaunayGraph_<T, 2> &delaunay) {
        auto &junction_points = graph.junctions;
        auto &adjacent_sites = graph.adjacent_sites;
        auto &adjacent_ridges = graph.adjacent_ridges;

        //const auto &facets = delaunay.facets;

        /*if (graph.d[site].vertices.empty()) {
            size_t cur_side = side;
            do {
                const auto &face = facets[triangle];
                size_t neighbour_id = face.neighbour_facets[side];
                size_t neighbour_side = face.neighbour_sides[side];
                const auto &neighbour = facets[neighbour_id];
                size_t neighbour_site = neighbour.vertices[neighbour_side];

                regions[site].junctions.push_back(triangle);
                regions[site].neighbour_regions.push_back(neighbour_site);

                size_t next_side = vx::NEXT_PTR_3[face.neighbour_sides[vx::NEXT_PTR_3[cur_side]]];
                triangle = face.neighbour_facets[vx::NEXT_PTR_3[cur_side]];
                cur_side = next_side;
            } while (triangle != vx::NO_PTR && regions[site].junctions[0] != triangle);

            if (triangle == vx::NO_PTR)
                regions[site].vertices = {vx::NO_PTR};
        }*/
        return true;
    }

    template<typename T>
    const VoronoiGraph<T, 2> &Voronoi2D<T>::getGraph() const {
        return graph;
    }

} // namespace vx

#endif //VORONOIX_VORONOID_2D_HPP
