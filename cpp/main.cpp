#include <iostream>
#include <vx/geometry/vector.hpp>
#include <vx/geometry/matrix.hpp>
#include <vx/geometry/polyhedron.hpp>
#include <vector>
#include <vx/algo/alg_delaunay_2d.hpp>
#include <xmmintrin.h>
#include <vx/algo/alg_voronoi_2d.hpp>
#include "vx/algo/alg_convex_hull_kd.hpp"


template<typename T, size_t DIM>
std::vector<vx::Vec<T, DIM>> select(std::vector<vx::Vec<T, DIM>> vs, size_t target_size) {
    std::vector<vx::Vec<T, DIM>> vsr;
    for (int i = 0; i < target_size; i++) {
        std::vector<vx::Vec<T, DIM>> candidates;
        for (int j = 0; j < 16; j++) {
            candidates.push_back(vs[rand() % vs.size()]);
        }
        vx::Vec<T, DIM> best;
        vx::Double best_dist = 0.0;
        for (auto c : candidates) {
            vx::Double nearest = 100000000.0;
            for (auto v : vsr) {
                nearest = std::min(nearest, vx::length(v - c));
            }
            if (nearest > best_dist) {
                best_dist = nearest;
                best = c;
            }
        }
        vsr.push_back(best);
    }
    return vsr;
}

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
    add(500, 0, 0);
    add(120, 15, -300);
    add(-300, 25, 180);
    add(250, -190, 29);*/
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
            auto offset = dot(mid, normal) - 8;
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

void cube_delaunay() {
    using T = long double;
    const size_t DIM = 3;

    std::vector<vx::Vec<T, 3>> datar;
    auto add = [&](T x, T y, T z) {
        datar.emplace_back();
        datar.back()[0] = x;
        datar.back()[1] = y;
        datar.back()[2] = z;
    };

    for (size_t i = 0; i < 128; i++) {
        add(rand() % 1900 - 950, rand() % 1900 - 950, rand() % 1900 - 950);
    }

    datar = select(datar, datar.size() / 4);

    std::vector<vx::Vec<T, 4>> data;
    for (auto &v: datar) {
        data.emplace_back();
        data.back()[0] = v[0];
        data.back()[1] = v[1];
        data.back()[2] = v[2];
        data.back()[3] = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / 3000;
    }

    vx::QuickHull<T, 4> voronoi_qhull(data);

    std::vector<std::vector<size_t>> neighbours(data.size());
    std::vector<vx::Polyhedron<T, 3>> polyhedrons;
    for (const auto &f : voronoi_qhull.getFacets()) {
        std::vector<vx::Vec<T, 3>> vecs;
        for (int i : f) vecs.push_back(datar[i]);
        polyhedrons.push_back(vx::makePolyhedron(vecs));
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
    auto &graph = delaunay.getGraph();
    for (size_t i = 0; i < graph.numFacets(); i++) {
        auto f = graph.getFacet(i);
        auto v = f.getVertices();
        out << "[" << v[0].getIndex() << ", " << v[1].getIndex() << ", " << v[2].getIndex() << "],\n";
    }
    out << "]}\n";
    out << "];\n";
}

void dump_convex(const std::vector<vx::Vec<vx::Double, 3>> &vs) {
    vx::QuickHull<vx::Double, 3> qhull(vs);
    std::ofstream out(R"(E:\dev\js\qhull-js\src\data.js)");
    out << "let polyhedrons = [\n";
    out << "{ \"v\":[\n";
    for (const auto &p : vs) {
        out << "[" << p[0] << ", " << p[1] << ", " << p[2] << "],\n";
    }
    out << "],\n";
    out << "\"f\":[\n";
    for (auto f : qhull.getFacets()) {
        out << "[" << f[0] << ", " << f[1] << ", " << f[2] << "],\n";
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
    /*vx::Voronoi2D<vx::Double> vor(vs);
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
    out << "];\n";*/
}

void make1() {


    std::vector<vx::Vec<vx::Double, 2>> vs;
    for (int i = 0; i < 32; i++) {
        std::vector<vx::Vec<vx::Double, 2>> candidates;
        for (int i = 0; i < 10; i++) {
            vx::Vec<vx::Double, 2> v(rand() % 4000 - 2000, rand() % 4000 - 2000);
            if (length(v) < 2000) {
                candidates.push_back(v);
            }
        }
        vx::Vec<vx::Double, 2> best;
        vx::Double best_dist = 0.0;
        for (auto c : candidates) {
            vx::Double nearest = 100000000.0;
            for (auto v : vs) {
                nearest = std::min(nearest, vx::length(v - c));
            }
            if (nearest > best_dist) {
                best_dist = nearest;
                best = c;
            }
        }
        vs.push_back(best);
    }

    dump_delaunay(vs);
}

void make2() {

    std::vector<vx::Vec<vx::Double, 3>> vs;
    for (int i = 0; i < 512; i++) {
        std::vector<vx::Vec<vx::Double, 3>> candidates;
        for (int i = 0; i < 10; i++) {
            vx::Vec<vx::Double, 3> v(rand() % 4000 - 2000, rand() % 4000 - 2000, rand() % 4000 - 2000);
            if (length(v) < 2000) {
                candidates.push_back(v);
            }
        }
        vx::Vec<vx::Double, 3> best;
        vx::Double best_dist = 0.0;
        for (auto c : candidates) {
            vx::Double nearest = 100000000.0;
            for (auto v : vs) {
                nearest = std::min(nearest, vx::length(v - c));
            }
            if (nearest > best_dist) {
                best_dist = nearest;
                best = c;
            }
        }
        vs.push_back(best);
    }

    dump_convex(vs);

}

void make_skip() {
    std::vector<vx::Vec<vx::Double, 2>> vs0;
    for (int i = 0; i < 32; i++) {
        std::vector<vx::Vec<vx::Double, 2>> candidates;
        for (int i = 0; i < 10; i++) {
            vx::Vec<vx::Double, 2> v(rand() % 4000 - 2000, rand() % 4000 - 2000);
            if (length(v) < 2000) {
                candidates.push_back(v);
            }
        }
        vx::Vec<vx::Double, 2> best;
        vx::Double best_dist = 0.0;
        for (auto c : candidates) {
            vx::Double nearest = 100000000.0;
            for (auto v : vs0) {
                nearest = std::min(nearest, vx::length(v - c));
            }
            if (nearest > best_dist) {
                best_dist = nearest;
                best = c;
            }
        }
        vs0.push_back(best);
    }
    std::vector<vx::Vec<vx::Double, 2>> vs1 = select(vs0, vs0.size() / 2);
    std::vector<vx::Vec<vx::Double, 2>> vs2 = select(vs1, vs1.size() / 2);
    std::vector<vx::Vec<vx::Double, 2>> vs3 = select(vs2, vs2.size() / 2);
    std::ofstream out(R"(E:\dev\js\qhull-js\src\data.js)");
    out << "let polyhedrons = [\n";
    for (auto vs : {vs0, vs1, vs2, vs3}) {
        vx::Delaunay2D<vx::Double> delaunay(vs);
        out << "{ \"v\":[\n";
        for (const auto &p : vs) {
            out << "[" << p[0] << ", " << p[1] << ", " << length_square(p) / 400000 << "],\n";
        }
        out << "],\n";
        out << "\"f\":[\n";
        auto &graph = delaunay.getGraph();
        for (size_t i = 0; i < graph.numFacets(); i++) {
            auto f = graph.getFacet(i);
            auto v = f.getVertices();
            out << "[" << v[0].getIndex() << ", " << v[1].getIndex() << ", " << v[2].getIndex() << "],\n";
        }
        out << "]},\n";
    }
    out << "];\n";
}

void make_convex() {

    std::vector<vx::Vec<vx::Double, 3>> vs;
    for (int i = 0; i < 512; i++) {
        std::vector<vx::Vec<vx::Double, 3>> candidates;
        for (int i = 0; i < 10; i++) {
            vx::Vec<vx::Double, 3> v(rand() % 4000 - 2000, rand() % 4000 - 2000, rand() % 4000 - 2000);
            if (length(v) < 2000) {
                candidates.push_back(v);
            }
        }
        vx::Vec<vx::Double, 3> best;
        vx::Double best_dist = 0.0;
        for (auto c : candidates) {
            vx::Double nearest = 100000000.0;
            for (auto v : vs) {
                nearest = std::min(nearest, vx::length(v - c));
            }
            if (nearest > best_dist) {
                best_dist = nearest;
                best = c;
            }
        }
        vs.push_back(best);
    }

    dump_convex(vs);
}

#include <chrono>

void bench2(int t = 0) {
    using namespace std::chrono;
    std::vector<vx::Vec<double, 2>> points;
    int n = 1'000'000;
    //srand(1);
    for(int i = 0; i < n; i++) {
        if (!t) {
            double x = rand() * 1.0 / RAND_MAX * 2000 - 1000;
            double y = rand() * 1.0 / RAND_MAX * 2000 - 1000;
            points.emplace_back(x, y);
        } else {
            double u = (rand() * 1.0 * RAND_MAX + rand()) * 1.0 / RAND_MAX / RAND_MAX * 0.88 + 0.1;
            double v = (rand() * 1.0 * RAND_MAX + rand()) * 1.0 / RAND_MAX / RAND_MAX;

            double x = sqrt(-2 * log(u)) * cos(2 * 3.141592 * v) * 1000;
            double y = sqrt(-2 * log(u)) * sin(2 * 3.141592 * v) * 1000;
            points.emplace_back(x, y);
        }
    }


    high_resolution_clock::time_point start = high_resolution_clock::now();
    for(int i = 0; i < 10; i++) {
        //vx::Voronoi2D<double> del(points);
        vx::GiftWrapping<double> ch(points);
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);

    std::cout << "It took me " << time_span.count() << " seconds.";
    std::cout << std::endl;
}

void bench3(int t = 0) {
    using namespace std::chrono;
    std::vector<vx::Vec<double, 3>> points;
    int n = 10'000;
    srand(0);
    for(int i = 0; i < n; i++) {
        if (!t) {
            double x = rand() * 1.0 / RAND_MAX * 2000 - 1000;
            double y = rand() * 1.0 / RAND_MAX * 2000 - 1000;
            double z = rand() * 1.0 / RAND_MAX * 2000 - 1000;
            points.emplace_back();
            points.back()[0] = x;
            points.back()[1] = y;
            points.back()[2] = z;
        } else {
            points.emplace_back();
            for(int j = 0; j < 3; j++) {
                double u = (rand() * 1.0 * RAND_MAX + rand()) * 1.0 / RAND_MAX / RAND_MAX * 0.88 + 0.1;
                double v = (rand() * 1.0 * RAND_MAX + rand()) * 1.0 / RAND_MAX / RAND_MAX;

                points.back()[j] = sqrt(-2 * log(u)) * cos(2 * 3.141592 * v) * 1000;
            }
        }
    }


    high_resolution_clock::time_point start = high_resolution_clock::now();
    for(int i = 0; i < 20; i++) {
        //vx::Voronoi2D<double> del(points);
        /*std::vector<vx::Vec<double, 4>> pts;
        pts.reserve(points.size());
        for(auto v : points) {
            vx::Vec<double, 4> u;
            u[0] = v[0];
            u[1] = v[1];
            u[2] = v[2];
            u[3] = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) / 3000;
            pts.push_back(u);
        }
        vx::QuickHull<double, 4> qhull(pts);*/
        vx::QuickHull<double, 3> qhull(points);
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);

    std::cout << "It took me " << time_span.count()/20 << " seconds.";
    std::cout << std::endl;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    cube_voronoi();
    std::cout << "Bye, World!" << std::endl;
}