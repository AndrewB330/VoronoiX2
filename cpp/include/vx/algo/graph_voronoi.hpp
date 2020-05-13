#ifndef VORONOIX_GRAPH_VORONOI_HPP
#define VORONOIX_GRAPH_VORONOI_HPP

#include <vx/geometry/vector.hpp>
#include <vector>
#include <array>

namespace vx {

    template<typename T>
    class Voronoi2D;

    template<typename T, size_t DIM>
    class VoronoiSite;

    template<typename T, size_t DIM>
    class VoronoiRidge;

    template<typename T, size_t DIM>
    class VoronoiGraph;

    template<typename T, size_t DIM>
    class VoronoiJunction;

    template<typename T, size_t DIM>
    class VoronoiGraph {
    public:
        VoronoiGraph() = default;

        VoronoiSite<T, DIM> getSite(size_t index) const;

        VoronoiJunction<T, DIM> getJunction(size_t index) const;

        size_t numSites() const;

        size_t numJunctions() const;

        friend vx::VoronoiSite<T, DIM>;
        friend vx::VoronoiRidge<T, DIM>;
        friend vx::VoronoiJunction<T, DIM>;

        friend vx::Voronoi2D<T>;

    protected:
        std::vector<vx::Vec<T, DIM>> points;
        std::vector<vx::Vec<T, DIM>> junctions;

        std::vector<std::vector<size_t>> ridges;
        std::vector<uint8_t> is_infinite; // bools

        std::vector<std::vector<size_t>> adjacent_sites;
        std::vector<std::vector<size_t>> adjacent_ridges;
    };

    template<typename T, size_t DIM>
    class VoronoiSite {
    public:
        bool isInfinite() const;

        size_t getIndex() const;

        vx::Vec<T, DIM> getPoint() const;

        std::vector<VoronoiSite> getAdjacent() const;

        std::vector<vx::VoronoiRidge<T, DIM>> getAdjacentRidges() const;

        friend vx::VoronoiGraph<T, DIM>;

    protected:
        VoronoiSite(size_t site_index, const vx::VoronoiGraph<T, DIM> &g) : g(g), site_index(site_index) {}

        const vx::VoronoiGraph<T, DIM> &g;
        size_t site_index;
    };

    template<typename T, size_t DIM>
    class VoronoiRidge {
    public:
        bool isInfinite() const;

        std::vector<vx::Vec<T, DIM>> getPoints() const;

        std::vector<vx::VoronoiJunction<T, DIM>> getJunctions() const;

        vx::VoronoiSite<T, DIM> getCurrentSite() const;

        vx::VoronoiSite<T, DIM> getNextSite() const;

        friend vx::VoronoiSite<T, DIM>;
        friend vx::VoronoiJunction<T, DIM>;

        friend vx::VoronoiGraph<T, DIM>;

    protected:
        VoronoiRidge(size_t site_index, size_t adjacent_index, const vx::VoronoiGraph<T, DIM> &g)
                : g(g), curr_site(site_index),
                  next_site(g.adjacent_sites[site_index][adjacent_index]),
                  ridge_index(g.adjacent_ridges[site_index][adjacent_index]) {
            VX_ASSERT(curr_site != next_site, "Sites should be different");
        }

        bool isReversed();

        const vx::VoronoiGraph<T, DIM> &g;
        size_t ridge_index;
        size_t curr_site;
        size_t next_site;
    };

    template<typename T, size_t DIM>
    class VoronoiJunction {
    public:
        size_t getIndex() const;

        vx::Vec<T, DIM> getPoint() const;

        friend vx::VoronoiSite<T, DIM>;
        friend vx::VoronoiRidge<T, DIM>;

        friend vx::VoronoiGraph<T, DIM>;
    protected:
        VoronoiJunction(size_t junction_index, const vx::VoronoiGraph<T, DIM> &g) : junction_index(junction_index),
                                                                                    g(g) {}

        size_t junction_index;
        const vx::VoronoiGraph<T, DIM> &g;
    };

    template<typename T, size_t DIM>
    size_t VoronoiJunction<T, DIM>::getIndex() const {
        return junction_index;
    }

    template<typename T, size_t DIM>
    Vec <T, DIM> VoronoiJunction<T, DIM>::getPoint() const {
        return g.junctions[junction_index];
    }

    /*
     * ========================================================
     * ==================== IMPLEMENTATION ====================
     * ========================================================
     */

    template<typename T, size_t DIM>
    bool VoronoiRidge<T, DIM>::isInfinite() const {
        return g.is_infinite[ridge_index];
    }

    template<typename T, size_t DIM>
    bool VoronoiRidge<T, DIM>::isReversed() {
        return curr_site > next_site;
    }

    template<typename T, size_t DIM>
    VoronoiSite<T, DIM> VoronoiRidge<T, DIM>::getNextSite() const {
        return VoronoiSite<T, DIM>(next_site, g);
    }

    template<typename T, size_t DIM>
    VoronoiSite<T, DIM> VoronoiRidge<T, DIM>::getCurrentSite() const {
        return VoronoiSite<T, DIM>(curr_site, g);
    }

    template<typename T, size_t DIM>
    std::vector<vx::Vec<T, DIM>> VoronoiRidge<T, DIM>::getPoints() const {
        std::vector<Vec<T, DIM>> res;
        res.reserve(g.ridges[ridge_index].size());
        for (size_t junction_index : g.ridges[ridge_index]) {
            res.push_back(g.junctions[junction_index]);
        }
        return res;
    }

    template<typename T, size_t DIM>
    std::vector<VoronoiJunction<T, DIM>> VoronoiRidge<T, DIM>::getJunctions() const {
        std::vector<VoronoiJunction<T, DIM>> res;
        res.reserve(g.ridges[ridge_index].size());
        for (size_t junction_index : g.ridges[ridge_index]) {
            res.push_back(VoronoiJunction<T, DIM>(junction_index, g));
        }
        return res;
    }

    template<typename T, size_t DIM>
    std::vector<vx::VoronoiSite<T, DIM>> VoronoiSite<T, DIM>::getAdjacent() const {
        std::vector<VoronoiSite> adjacent_sites;
        adjacent_sites.reserve(g.adjacent_sites[site_index].size());

        for (size_t adj_site : g.adjacent_sites[site_index]) {
            adjacent_sites.push_back(VoronoiSite(adj_site, g));
        }

        return adjacent_sites;
    }

    template<typename T, size_t DIM>
    Vec <T, DIM> VoronoiSite<T, DIM>::getPoint() const {
        return g.points[site_index];
    }

    template<typename T, size_t DIM>
    size_t VoronoiSite<T, DIM>::getIndex() const {
        return site_index;
    }

    template<typename T, size_t DIM>
    std::vector<VoronoiRidge<T, DIM>> VoronoiSite<T, DIM>::getAdjacentRidges() const {
        std::vector<VoronoiRidge<T, DIM>> adjacent_ridges;
        adjacent_ridges.reserve(g.adjacent_ridges[site_index].size());

        for (size_t i = 0; i < g.adjacent_ridges[site_index].size(); i++) {
            adjacent_ridges.push_back(VoronoiRidge<T, DIM>(site_index, i, g));
        }

        return adjacent_ridges;
    }

    template<typename T, size_t DIM>
    bool VoronoiSite<T, DIM>::isInfinite() const {
        for (size_t adj_ridge : g.adjacent_ridges[site_index]) {
            if (g.is_infinite[adj_ridge]) {
                return true;
            }
        }
        return false;
    }

    template<typename T, size_t DIM>
    VoronoiSite<T, DIM> VoronoiGraph<T, DIM>::getSite(size_t index) const {
        return vx::VoronoiSite(index, *this);
    }

    template<typename T, size_t DIM>
    size_t VoronoiGraph<T, DIM>::numSites() const {
        return points.size();
    }

    template<typename T, size_t DIM>
    size_t VoronoiGraph<T, DIM>::numJunctions() const {
        return junctions.size();
    }

    template<typename T, size_t DIM>
    VoronoiJunction<T, DIM> VoronoiGraph<T, DIM>::getJunction(size_t index) const {
        return VoronoiJunction<T, DIM>(index, *this);
    }

} // namespace vx

#endif //VORONOIX_GRAPH_VORONOI_HPP
