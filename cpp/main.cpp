#include <iostream>
#include <vx/geometry/vector.hpp>
#include <vx/geometry/matrix.hpp>
#include <vx/geometry/polyhedron.hpp>
#include <vector>
#include <vx/algo/vx_builder_2d.hpp>
#include <xmmintrin.h>
#include <vx/algo/alg_voronoi_2d.hpp>
#include "vx/algo/alg_convex_hull_kd.hpp"

template<typename T, size_t DIM>
void dump(const std::vector<vx::Polyhedron<T, DIM>> &polyhedrons) {
    std::ofstream out(R"(E:\dev\js\qhull-js\src\data.js)");
    out << "let polyhedrons = [\n";
    size_t offset = 0;
    for (auto poly : polyhedrons) {
        out << "{ \"v\":[\n";
        for (auto p : poly.points) {
            out << "[" << p[0] << "," << p[1] << "," << p[2] << "]," << std::endl;
        }
        out << "], \"f\": [\n";
        for (auto f : poly.facets) {
            //if (f.half_space.normal[DIM - 1] < 0)
            out << "[" << f[0] + offset << ", " << f[1] + offset << ", " << f[2] + offset
                << "],"
                << std::endl;
        }
        out << "]},";
        //offset += poly.points.size();
    }
    out << "];\n";
}

void cube_voronoi() {
    using T = long double;
    const size_t DIM = 3;

    std::vector<vx::Vec<T, 4>> data;
    std::vector<vx::Vec<T, 3>> datar;
    auto add = [&](T x, T y, T z) {
        data.emplace_back();
        data.back()[0] = x;
        data.back()[1] = y;
        data.back()[2] = z;
        data.back()[3] = (x * x + y * y + z * z) / 3000;
    };
    /*add(-500, 0, 0);
    add(500, 0, 0);*/
    /*auto gen = []() {
        int a = rand() % 1900 - 950;
        int b = rand() % 1900 - 950;
        if (abs(a) > abs(b)) return a;
        return b;
    };*/
    for (size_t i = 0; i < 150; i++) {
        add(rand() % 1900 - 950, rand() % 1900 - 950, rand() % 1900 - 950);
    }

    std::sort(data.begin(), data.end());
    for (auto v : data) datar.emplace_back(v[0], v[1], v[2]);

    vx::QuickHull<T, 4> voronoi_qhull(data);

    std::vector<std::vector<size_t>> neighbours(data.size());
    for (const auto &f : voronoi_qhull.getFacets()) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                if (i != j)
                    neighbours[f[i]].push_back(f[j]);
            }
        }
    }
    for (size_t i = 0; i < data.size(); i++) {
        std::sort(neighbours[i].begin(), neighbours[i].end());
        auto it = std::unique(neighbours[i].begin(), neighbours[i].end());
        neighbours[i].erase(it, neighbours[i].end());
    }

    std::vector<vx::Polyhedron<T, DIM>> polyhedrons;
    for (size_t i = 0; i < data.size(); i++) {
        std::vector<vx::HalfSpace<T, DIM>> hs;
        for (size_t j = 0; j < data.size(); j++) {
            //for (size_t j : neighbour_facets[i]) {
            //VX_ASSERT(i != j, "er");
            if (i == j) continue;
            auto mid = (datar[j] + datar[i]) * 0.5;
            auto normal = datar[j] - datar[i];
            normal /= length(normal);
            auto offset = dot(mid, normal) - 4;
            hs.emplace_back(normal, offset);
        }
        hs.emplace_back(vx::Vec<T, DIM>(1, 0, 0), 1000);
        hs.emplace_back(vx::Vec<T, DIM>(-1, 0, 0), 1000);
        hs.emplace_back(vx::Vec<T, DIM>(0, 1, 0), 1000);
        hs.emplace_back(vx::Vec<T, DIM>(0, -1, 0), 1000);
        hs.emplace_back(vx::Vec<T, DIM>(0, 0, 1), 1000);
        hs.emplace_back(vx::Vec<T, DIM>(0, 0, -1), 1000);
        polyhedrons.push_back(vx::makePolyhedron(hs, datar[i]));
    }
    dump(polyhedrons);
}

void sphere_voronoi() {
    using T = long double;
    const size_t DIM = 3;

    std::vector<vx::Vec<T, 4>> data;
    auto add = [&](T x, T y, T z) {
        if (sqrt(x * x + y * y + z * z) > 950) return;
        data.emplace_back();
        data.back()[0] = x;
        data.back()[1] = y;
        data.back()[2] = z;
        data.back()[3] = (x * x + y * y + z * z) / 3000;
    };
    for (size_t i = 0; i < 900; i++) {
        add(rand() % 1900 - 950, rand() % 1900 - 950, rand() % 1900 - 950);
    }
    std::vector<vx::Vec<T, 3>> datar(data.size());

    std::sort(data.begin(), data.end());
    size_t IT = 500;
    for (size_t k = 0; k < IT; k++) {
        for (size_t i = 0; i < data.size(); i++) datar[i] = vx::Vec<T, DIM>(data[i][0], data[i][1], data[i][2]);

        vx::QuickHull<T, 4> voronoi_qhull(data);

        std::vector<std::vector<size_t>> neighbours(data.size());
        for (const auto &f : voronoi_qhull.getFacets()) {
            for (size_t i = 0; i < 4; i++) {
                for (size_t j = 0; j < 4; j++) {
                    if (i != j)
                        neighbours[f[i]].push_back(f[j]);
                }
            }
        }
        for (size_t i = 0; i < data.size(); i++) {
            std::sort(neighbours[i].begin(), neighbours[i].end());
            auto it = std::unique(neighbours[i].begin(), neighbours[i].end());
            neighbours[i].erase(it, neighbours[i].end());
        }
        std::vector<vx::HalfSpace<T, DIM>> rhs;
        for (size_t i = 0; i < 0; i++) {
            int a = rand();
            int b = rand();
            double x = sin(a);
            double y = cos(a) * sin(b);
            double z = cos(a) * cos(b);
            rhs.emplace_back(vx::Vec<T, DIM>(x, y, z), 1000);
        }
        rhs.emplace_back(vx::Vec<T, DIM>(1, 0, 0), 1000);
        rhs.emplace_back(vx::Vec<T, DIM>(-1, 0, 0), 1000);
        rhs.emplace_back(vx::Vec<T, DIM>(0, 1, 0), 1000);
        rhs.emplace_back(vx::Vec<T, DIM>(0, -1, 0), 1000);
        rhs.emplace_back(vx::Vec<T, DIM>(0, 0, 1), 1000);
        rhs.emplace_back(vx::Vec<T, DIM>(0, 0, -1), 1000);

        T s3 = std::sqrt(3);
        rhs.emplace_back(vx::Vec<T, DIM>(1, 1, 1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(1, 1, -1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(1, -1, 1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(1, -1, -1) / s3, 1400);

        rhs.emplace_back(vx::Vec<T, DIM>(-1, 1, 1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(-1, 1, -1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(-1, -1, 1) / s3, 1400);
        rhs.emplace_back(vx::Vec<T, DIM>(-1, -1, -1) / s3, 1400);

        std::vector<vx::Polyhedron<T, DIM>> polyhedrons;
        //polyhedrons.push_back(vx::makePolyhedron(rhs, vx::Vec<T, DIM>(0, 0, 0)));
        auto ndata = data;
        for (size_t i = 0; i < data.size(); i++) {
            std::vector<vx::HalfSpace<T, DIM>> hs;
            //for (size_t j = 0; j < data.size(); j++) {
            T sum = T(0);
            for (size_t j : neighbours[i]) {
                //VX_ASSERT(i != j, "er");
                if (i == j) continue;
                auto mid = (datar[j] + datar[i]) * 0.5;
                auto normal = datar[j] - datar[i];
                normal /= length(normal);
                auto offset = dot(mid, normal) - 4;
                hs.emplace_back(normal, offset);
            }
            hs.insert(hs.end(), rhs.begin(), rhs.end());
            auto p = vx::makePolyhedron(hs, datar[i]);
            polyhedrons.push_back(p);
            vx::Vec<T, DIM> c;
            for (auto f : p.facets) {
                vx::Mat<T, DIM, DIM> m;
                m.set_row(0, p.points[f[0]]);
                m.set_row(1, p.points[f[1]]);
                m.set_row(2, p.points[f[2]]);
                T area = std::abs(det(m));
                for (size_t j = 0; j < DIM; j++) {
                    c += p.points[f[j]] * area;
                    sum += area;
                }
            }
            c /= sum;
            c = (c + c) * 0.5;
            vx::Vec<T, 4> c4;
            c4[0] = c[0];
            c4[1] = c[1];
            c4[2] = c[2];
            c4[3] = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / 3000;
            ndata[i] = c4;
        }
        data = ndata;
        if (k % 10 == 0)
            dump(polyhedrons);
        std::cout << k << std::endl;
    }
}

void dump_delaunay(const std::vector<vx::Vec<vx::Double, 2>> &vs) {
    vx::Delaunay2D<vx::Double> delaunay(vs);
    std::ofstream out(R"(E:\dev\js\qhull-js\src\data.js)");
    out << "let polyhedrons = [\n";
    out << "{ \"v\":[\n";
    for (const auto &p : vs) {
        out << "[" << p[0] << ", " << p[1] << ", " << length_square(p) / 400000 << "],\n";
    }
    out << "],\n";
    out << "\"f\":[\n";
    auto & graph = delaunay.getGraph();
    for (size_t i = 0; i < graph.numFacets(); i++) {
        auto f = graph.getFacet(i);
        auto v = f.getVertices();
        out << "[" << v[0].getIndex() << ", " << v[1].getIndex() << ", " << v[2].getIndex() << "],\n";
    }
    out << "]}\n";
    out << "];\n";
}

void dump_voronoi(std::vector<vx::Vec<vx::Double, 2>> vs) {
    /*for (size_t k = 0; k < 80; k++) {
        std::cout << k << std::endl;
        vx::Voronoi2D<vx::Double> vor(vs);
        auto graph = vor.getGraph();
        auto nvs = std::vector<vx::Vec<vx::Double, 2>>(vs.size());
        for (size_t i = 0; i < vs.size(); i++) {
            bool ok = !vor.regions[i].vertices.empty();
            for (size_t c : vor.regions[i].vertices) {
                if (std::abs(vor.centers[c].x) > 5000) ok = false;
                if (std::abs(vor.centers[c].y) > 5000) ok = false;
                nvs[i] += vor.centers[c];
            }
            if (!ok)
                nvs[i] = vs[i];
            else
                nvs[i] /= vor.regions[i].vertices.size();
        }
        vs = nvs;
    }*/
    vx::Voronoi2D<vx::Double> vor(vs);
    const auto & graph = vor.getGraph();
    std::ofstream out(R"(E:\dev\js\qhull-js\src\data.js)");
    out << "let polyhedrons = [\n";
    out << "{ \"v\":[\n";
    for (size_t i = 0; i < graph.numJunctions(); i++) {
        auto p = graph.getJunction(i).getPoint();
        out << "[" << p[0] << ", " << p[1] << ", " << -length_square(p) / 6000 << "],\n";
    }
    for (const auto &p : vs) {
        out << "[" << p[0] << ", " << p[1] << ", " << -length_square(p) / 6000 << "],\n";
    }
    out << "],\n";
    out << "\"f\":[\n";
    for (size_t i = 0; i < vs.size(); i++) {
        auto v = graph.getSite(i);

        for (auto r : v.getAdjacentRidges()) {
            auto pts = r.getJunctions();
            out << "[" << pts[0].getIndex() << ", " << pts[1].getIndex() << ", " << i + graph.numSites() << "],\n";
        }
    }
    out << "],\n";
    out << "\"c\":[\n";
    for (size_t i = 0; i < vs.size(); i++) {
        size_t ra = rand() % 16;
        size_t rb = rand() % 16;
        size_t ga = rand() % 16;
        size_t gb = rand() % 16;
        size_t ba = rand() % 16;
        size_t bb = rand() % 16;
        //out << "0x" << std::hex << ra << rb << ga << gb << ba << bb << ",\n";
        size_t adjacent_num = graph.getSite(i).getAdjacent().size();
        for (size_t j = 0; j < adjacent_num; j++) {
            out << "0x" << std::hex << ra << rb << ga << gb << ba << bb << ",\n";
        }
    }
    out << std::dec;
    out << "]}\n";
    out << "];\n";
}

int main() {
    std::cout << "Hello, World!" << std::endl;

    std::vector<vx::Vec<vx::Double, 2>> vs;
    for(int i = 0; i < 1000; i++) {
        vs.emplace_back(rand()%2000-1000, rand()%2000-1000);
    }

    dump_voronoi(vs);

    std::cout << "Bye, World!" << std::endl;
}