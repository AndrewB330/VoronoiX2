#include <iostream>
#include <vx/math/vec.hpp>
#include <vx/math/mat.hpp>
#include <vx/math/qhull.hpp>
#include <vector>
#include <vx/math/geometry.hpp>

void test_2() {
    std::vector<vx::Vec<double, 2>> v;
    //v.emplace_back(0, 0);
    v.emplace_back(3, 3);
    v.emplace_back(-3, 2);
    v.emplace_back(1, -3);
    //v.emplace_back(1.1, -2);
    v.emplace_back(1.5, -2.5);
    v.emplace_back(0.5, 3.2);
    //vx::Context pts(v);
    vx::QHull q_hull(v);
}

void test_3() {
    std::vector<vx::Vec<long double, 3>> v;
    v.emplace_back(1, 1, 1);
    v.emplace_back(-1, 1, 1);
    //v.emplace_back(0, 0, 0);
    v.emplace_back(0, -1, 2);
    v.emplace_back(0, -1, -2);
    v.emplace_back(-0.52, 0, 0);
    v.emplace_back(0.7, -0.2, -0.1);
    v.emplace_back(0, 0.6, -0.2);
    //vx::Context pts(v);
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
}

void test_3_r() {
    std::vector<vx::Vec<long double, 3>> v;
    for (size_t i = 0; i < 128; i++) {
        double a = rand() % 331;
        double b = rand() % 256;
        double y = sin(a) * 1024 + (rand() % 100 / 100.0);
        double x = cos(b) * cos(a) * 1024 + (rand() % 100 / 100.0);
        double z = sin(b) * cos(a) * 1024 + (rand() % 100 / 100.0);
        v.emplace_back(x, y, z);
        //std::cout << v.back() << std::endl;
    }
    /*for (size_t i = 0; i < 100; i++) {
        vx::QHull q_hull(v);
        std::vector<vx::Vec<long double, 3>> acc(v.size());
        std::vector<int> cnt(v.size());
        for (const auto &f : q_hull.facets) {
            for (size_t j = 0; j < 3; j++) {
                size_t a = f->vertices[j];
                size_t b = f->vertices[(j + 1) % 3];
                acc[a] += v[b];
                cnt[a]++;
            }
        }
        for (size_t j = 0; j < v.size(); j++) {
            v[j] = v[j] * 0.1 + acc[j] / cnt[j] * 0.9;
        }
        for (auto &p : v) {
            p = p * (1024 / length(p));
        }
    }*/
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json("E:\\dev\\js\\qhull-js\\src\\data.js");
}

void test_3_v() {
    std::set<vx::Vec<long double, 3>> v;
    for (size_t i = 0; i < 1024 * 128; i++) {
        double x = rand() % 2001 - 1000;
        double y = rand() % 2001 - 1000;
        vx::Vec<long double, 3> u;
        u[0] = x;
        u[1] = y;
        u[2] = (x * x + y * y) / 2000;
        v.insert(u);
    }
    std::cout << v.size() << std::endl;
    std::vector<vx::Vec<long double, 3>> vv(v.begin(), v.end());
    vx::QHull q_hull(vv);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json(R"(E:\dev\js\qhull-js\src\data.js)");
}

void test_4_v() {
    std::set<vx::Vec<long double, 4>> v;
    for (size_t i = 0; i < 1024 * 32; i++) {
        double x = rand() % 2001 - 1000;
        double y = rand() % 2001 - 1000;
        double z = rand() % 2001 - 1000;
        vx::Vec<long double, 4> u;
        u[0] = x;
        u[1] = y;
        u[2] = z;
        u[3] = (x * x + y * y + z * z) / 3000;
        v.insert(u);
    }
    std::cout << v.size() << std::endl;
    std::vector<vx::Vec<long double, 4>> vv(v.begin(), v.end());
    vx::QHull q_hull(vv);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json(R"(E:\dev\js\qhull-js\src\data.js)");
}

void test_2_r() {
    std::vector<vx::Vec<long double, 2>> v;
    size_t n = 2000;
    for (size_t i = 0; i < n; i++) {
        double a = i * 3.14 * 2 / n;
        double dx = (rand() % 100) / 10000.0;
        double dy = (rand() % 100) / 10000.0;
        v.emplace_back(cos(a) * 1024, sin(a) * 1024);
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
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
            out << "[" << f.vertices[0] + offset << ", " << f.vertices[1] + offset << ", " << f.vertices[2] + offset
                << "],"
                << std::endl;
        }
        out << "]},";
        //offset += poly.points.size();
    }
    out << "];\n";
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    /*std::vector<vx::Vec<double, 2>> v;
    v.emplace_back(3, 0);
    v.emplace_back(3, 3);
    v.emplace_back(0, 5);
    v.emplace_back(-4, 2);
    v.emplace_back(-3, -3);
    v.emplace_back(1.5, -1.5);
    for (auto &v : v) {
        v /= length_square(v);
    }
    vx::QHull<double, 2> hull(v);
    for (auto f : hull.facets) {
        auto v = f->half_space.normal * f->half_space.offset;
        std::cout << v / length_square(v) << std::endl;
    }*/
    using T = long double;
    const size_t DIM = 3;
    vx::GeometryContext<T> ctx(1000.0);

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
    for(size_t i = 0; i < 195; i++) {
        add(rand() % 1900 - 950, rand() % 1900 - 950, rand() % 1900 - 950);
    }

    std::sort(data.begin(), data.end());
    for(auto v : data) datar.emplace_back(v[0], v[1], v[2]);

    vx::QHull<T, 4> voronoi_qhull(data);

    std::vector<std::vector<size_t>> neighbours(data.size());
    for (auto f : voronoi_qhull.getHull().facets) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                if (i != j)
                    neighbours[f.vertices[i]].push_back(f.vertices[j]);
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
        for (size_t j : neighbours[i]) {
            VX_ASSERT(i != j, "er");
            if (i == j) continue;
            auto mid = (datar[j] + datar[i]) * 0.5;
            auto normal = datar[j] - datar[i];
            normal /= length(normal);
            auto offset = dot(mid, normal) - 25;
            hs.emplace_back(normal, offset, ctx);
        }
        hs.emplace_back(vx::Vec<T, DIM>(1, 0, 0), 1000, ctx);
        hs.emplace_back(vx::Vec<T, DIM>(-1, 0, 0), 1000, ctx);
        hs.emplace_back(vx::Vec<T, DIM>(0, 1, 0), 1000, ctx);
        hs.emplace_back(vx::Vec<T, DIM>(0, -1, 0), 1000, ctx);
        hs.emplace_back(vx::Vec<T, DIM>(0, 0, 1), 1000, ctx);
        hs.emplace_back(vx::Vec<T, DIM>(0, 0, -1), 1000, ctx);
        polyhedrons.push_back(vx::makePolyhedron(hs, datar[i]));
    }
    dump(polyhedrons);

    std::cout << "DONE" << std::endl;
}