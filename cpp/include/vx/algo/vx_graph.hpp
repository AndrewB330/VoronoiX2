#ifndef VORONOIX_VX_GRAPH_HPP
#define VORONOIX_VX_GRAPH_HPP

#include <vector>
#include <vx/geometry/vector.hpp>

namespace vx {

    template<typename T, size_t DIM>
    struct VxGraph {
        std::vector<Vec<T, DIM>> points;
        std::vector<Vec<T, DIM>> joints;

        std::vector<std::array<uint64_t, DIM + 1>> facet_vertices;
        std::vector<std::array<uint64_t, DIM + 1>> facet_adjacent;
        std::vector<std::array<uint8_t, DIM + 1>> facet_adjacent_side;

        std::vector<std::vector<uint64_t>> vertex_adjacent;
        std::vector<std::vector<uint64_t>> vertex_facets;
        std::vector<bool> vertex_boundary;

        std::vector<std::vector<std::vector<uint64_t>>> edge_joints; // well... it's ok

        uint64_t boundary_vertex;
    };

} // namespace vx

#endif //VORONOIX_VX_GRAPH_HPP
