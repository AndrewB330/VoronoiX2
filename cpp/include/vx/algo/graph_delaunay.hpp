#ifndef VORONOIX_GRAPH_DELAUNAY_HPP
#define VORONOIX_GRAPH_DELAUNAY_HPP

#include <vx/geometry/vector.hpp>
#include <vector>
#include <array>
#include "alg_voronoi_kd.hpp"

namespace vx {

    template<typename T>
    class Delaunay2D;

    template<typename T>
    class Voronoi2D;

    template<typename T, size_t DIM>
    class DelaunayVertex;

    template<typename T, size_t DIM>
    class DelaunayFacet_;

    template<typename T, size_t DIM>
    class DelaunayEdge;

    template<typename T, size_t DIM>
    class DelaunayGraph_;

    template<typename T, size_t DIM>
    class DelaunayGraph_ {
    public:
        DelaunayGraph_() = default;

        vx::DelaunayVertex<T, DIM> getVertex(size_t vertex) const;

        vx::DelaunayFacet_<T, DIM> getFacet(size_t facet) const;

        size_t numVertices() const;

        size_t numFacets() const;

        friend DelaunayVertex<T, DIM>;
        friend DelaunayFacet_<T, DIM>;
        friend DelaunayEdge<T, DIM>;

        friend Delaunay2D<T>;
        friend Voronoi2D<T>;
    protected:

        std::vector<vx::Vec<T, DIM>> points;
        std::vector<std::array<size_t, DIM + 1>> facets_vertices;
        std::vector<std::array<size_t, DIM + 1>> facets_adjacent;
        std::vector<std::array<size_t, DIM + 1>> facets_adjacent_side;

        std::vector<std::vector<size_t>> edge;

        std::vector<std::vector<size_t>> adjacent_facets;
        std::vector<std::vector<size_t>> adjacent_points;
        std::vector<std::vector<size_t>> adjacent_edges;
    };

    template<typename T, size_t DIM>
    class DelaunayVertex {
    public:
        size_t getIndex() const;

        const vx::Vec<T, DIM> &getPoint() const;

        std::vector<vx::DelaunayEdge<T, DIM>> getAdjacentEdges() const;

        std::vector<vx::DelaunayVertex<T, DIM>> getAdjacentVertices() const;

        std::vector<vx::DelaunayFacet_<T, DIM>> getAdjacentFacets() const;

        friend DelaunayGraph_<T, DIM>;
        friend DelaunayFacet_<T, DIM>;
        friend DelaunayEdge<T, DIM>;

        friend Delaunay2D<T>;
        friend Voronoi2D<T>;
    protected:
        DelaunayVertex(size_t vertex, const DelaunayGraph_<T, DIM> &g) : vertex(vertex), g(g) {}

        size_t vertex;
        const DelaunayGraph_<T, DIM> &g;
    };

    template<typename T, size_t DIM>
    class DelaunayEdge {
    public:
        vx::DelaunayVertex<T, DIM> getVertexFrom() const;

        vx::DelaunayVertex<T, DIM> getVertexTo() const;

        const vx::Vec<T, DIM> &getPointFrom() const;

        const vx::Vec<T, DIM> &getPointTo() const;

        //std::vector<vx::DelaunayFacet_<T, DIM>> getAdjacentFacets() const;

        friend DelaunayGraph_<T, DIM>;
        friend DelaunayFacet_<T, DIM>;
        friend DelaunayVertex<T, DIM>;

        friend Delaunay2D<T>;
    protected:
        DelaunayEdge(size_t vertex, size_t i, const DelaunayGraph_<T, DIM> &g) : vertex(vertex), adj_index(i), g(g) {}

        size_t vertex;
        size_t adj_index;
        const DelaunayGraph_<T, DIM> &g;
    };

    template<typename T, size_t DIM>
    class DelaunayFacet_ {
    public:

        vx::DelaunayVertex<T, DIM> getVertexOppositeToSide(size_t side) const;

        vx::DelaunayFacet_<T, DIM> getFacetAtSide(size_t side) const;

        std::vector<vx::DelaunayVertex<T, DIM>> getVertices() const;

        std::vector<vx::DelaunayFacet_<T, DIM>> getAdjacnetFacets() const;

        friend DelaunayGraph_<T, DIM>;
        friend DelaunayEdge<T, DIM>;
        friend DelaunayVertex<T, DIM>;

        friend Delaunay2D<T>;
    protected:
        DelaunayFacet_(size_t facet, size_t shift, const DelaunayGraph_<T, DIM> &g) : facet(facet), shift(shift),
                                                                                      g(g) {}

        size_t getShiftedSide(size_t side);

        size_t facet;
        size_t shift;
        const DelaunayGraph_<T, DIM> &g;
    };

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    DelaunayVertex<T, DIM> DelaunayGraph_<T, DIM>::getVertex(size_t vertex) const {
        return DelaunayVertex<T, DIM>(vertex, *this);
    }

    template<typename T, size_t DIM>
    size_t DelaunayGraph_<T, DIM>::numVertices() const {
        return points.size();
    }

    template<typename T, size_t DIM>
    size_t DelaunayGraph_<T, DIM>::numFacets() const {
        return facets_vertices.size();
    }

    template<typename T, size_t DIM>
    DelaunayFacet_<T, DIM> DelaunayGraph_<T, DIM>::getFacet(size_t facet) const {
        return DelaunayFacet_<T, DIM>(facet, 0, *this);
    }

    template<typename T, size_t DIM>
    const Vec <T, DIM> &DelaunayVertex<T, DIM>::getPoint() const {
        return g.points[vertex];
    }

    template<typename T, size_t DIM>
    std::vector<DelaunayEdge<T, DIM>> DelaunayVertex<T, DIM>::getAdjacentEdges() const {
        std::vector<DelaunayEdge<T, DIM>> adjacent;
        adjacent.reserve(g.adjacent_points[vertex].size());
        for (size_t i = 0; i < g.adjacent_points[vertex].size(); i++) {
            adjacent.push_back(DelaunayEdge<T, DIM>(vertex, i, g));
        }
        return adjacent;
    }

    template<typename T, size_t DIM>
    std::vector<DelaunayVertex<T, DIM>> DelaunayVertex<T, DIM>::getAdjacentVertices() const {
        std::vector<DelaunayVertex<T, DIM>> adjacent;
        adjacent.reserve(g.adjacent_points[vertex].size());
        for (size_t i = 0; i < g.adjacent_points[vertex].size(); i++) {
            adjacent.push_back(DelaunayVertex<T, DIM>(g.adjacent_points[i], g));
        }
        return adjacent;
    }

    template<typename T, size_t DIM>
    std::vector<DelaunayFacet_<T, DIM>> DelaunayVertex<T, DIM>::getAdjacentFacets() const {
        std::vector<DelaunayFacet_<T, DIM>> adjacent;
        adjacent.reserve(g.adjacent_facets[vertex].size());
        for (size_t i = 0; i < g.adjacent_facets[vertex].size(); i++) {
            size_t facet = g.adjacent_facets[i];
            for (size_t side = 0; side < (DIM + 1); side++) {
                if (g.facets_vertices[facet][side] == vertex) {
                    adjacent.push_back(DelaunayFacet_<T, DIM>(facet, side, g));
                }
            }
        }
        return adjacent;
    }

    template<typename T, size_t DIM>
    size_t DelaunayVertex<T, DIM>::getIndex() const {
        return vertex;
    }

    template<typename T, size_t DIM>
    DelaunayVertex<T, DIM> DelaunayEdge<T, DIM>::getVertexFrom() const {
        return DelaunayVertex<T, DIM>(vertex, g);
    }

    template<typename T, size_t DIM>
    DelaunayVertex<T, DIM> DelaunayEdge<T, DIM>::getVertexTo() const {
        return DelaunayVertex<T, DIM>(g.adjacent_points[vertex], g);
    }

    template<typename T, size_t DIM>
    const Vec <T, DIM> &DelaunayEdge<T, DIM>::getPointFrom() const {
        return g.points[vertex];
    }

    template<typename T, size_t DIM>
    const Vec <T, DIM> &DelaunayEdge<T, DIM>::getPointTo() const {
        return g.points[g.adjacent_points[vertex]];
    }

    template<typename T, size_t DIM>
    DelaunayVertex<T, DIM> DelaunayFacet_<T, DIM>::getVertexOppositeToSide(size_t side_shift) const {
        size_t actual_side = shift + side_shift;
        if (actual_side >= DIM + 1)
            actual_side -= DIM + 1;
        return DelaunayVertex<T, DIM>(g.facets_vertices[facet][actual_side], g);
    }

    template<typename T, size_t DIM>
    size_t DelaunayFacet_<T, DIM>::getShiftedSide(size_t side) {
        size_t res = side + shift;
        if (res >= DIM + 1)
            res -= DIM + 1;
        VX_ASSERT(res >= 0 && res < DIM + 1, "side out of range");
        return res;
    }

    template<typename T, size_t DIM>
    DelaunayFacet_<T, DIM> DelaunayFacet_<T, DIM>::getFacetAtSide(size_t side) const {
        return DelaunayFacet_<T, DIM>(g.facets_adjacent[getShiftedSide(side)],
                                      g.facets_adjacent_side[getShiftedSide(side)], g);
    }

    template<typename T, size_t DIM>
    std::vector<DelaunayVertex<T, DIM>> DelaunayFacet_<T, DIM>::getVertices() const {
        std::vector<DelaunayVertex<T, DIM>> vertices;
        for (size_t i = 0; i < DIM + 1; i++) {
            vertices.push_back(DelaunayVertex<T, DIM>(g.facets_vertices[facet][i], g));
        }
        return vertices;
    }

    template<typename T, size_t DIM>
    std::vector<DelaunayFacet_<T, DIM>> DelaunayFacet_<T, DIM>::getAdjacnetFacets() const {
        std::vector<DelaunayFacet_<T, DIM>> facets;
        for (size_t i = 0; i < DIM + 1; i++) {
            facets.push_back(getFacetAtSide(i));
        }
        return facets;
    }

} // namespace vx

#endif //VORONOIX_GRAPH_DELAUNAY_HPP
