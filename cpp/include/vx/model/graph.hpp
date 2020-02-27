#ifndef VORONOIX_GRAPH_HPP
#define VORONOIX_GRAPH_HPP

#include <vector>

namespace vx {

    template<typename E>
    struct Edge {
        int from;
        int to;
        E value;
    };

    template<typename V>
    struct Vertex {
        int id;
        V value;
    };

    template<typename V, typename E>
    class SimpleGraph {
    public:
        explicit SimpleGraph(const std::vector<V> &data, const std::vector<std::vector<std::pair<int, E>>> &edges);

        void add_edge(int from, int to, E data);

        Vertex<V> get_vertex(int index) const;

        const std::vector<Edge<E>> &get_edges(int index) const;

        unsigned int size() const;

    private:
        std::vector<std::vector<Edge<E>>> edges;
        std::vector<V> data;
    };

} // namespace vx

#endif //VORONOIX_GRAPH_HPP
