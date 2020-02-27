#include <iostream>
#include <vx/math/vec.hpp>
#include <vx/math/mat.hpp>
#include <vx/math/qhull.hpp>
#include <vector>

void test_2() {
    std::vector<vx::Vec<double, 2>> v;
    v.emplace_back(0, 0);
    v.emplace_back(3, 3);
    v.emplace_back(-3, 2);
    v.emplace_back(1, -3);
    v.emplace_back(1.1, -2);
    //vx::Context ctx(v);
    vx::QHull q_hull(v);
}

void test_3() {
    std::vector<vx::Vec<long double, 3>> v;
    v.emplace_back(1, 1, 1);
    v.emplace_back(-1, 1, 1);
    v.emplace_back(0, 0, 0);
    v.emplace_back(0, -1, 2);
    v.emplace_back(0, -1, -2);
    v.emplace_back(-0.52, 0, 0);
    //vx::Context ctx(v);
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
}

void test_3_r() {
    std::vector<vx::Vec<long double, 3>> v;
    for (size_t i = 0; i < 256; i++) {
        double a = rand() % 331;
        double b = rand() % 256;
        double y = sin(a) * 1024 + (rand() % 100 / 100.0);
        double x = cos(b) * cos(a) * 1024 + (rand() % 100 / 100.0);
        double z = sin(b) * cos(a) * 1024 + (rand() % 100 / 100.0);
        v.emplace_back(x, y, z);
        //std::cout << v.back() << std::endl;
    }
    for(size_t i = 0; i < 100; i++) {
        vx::QHull q_hull(v);
        std::vector<vx::Vec<long double, 3>> acc(v.size());
        std::vector<int> cnt(v.size());
        for(const auto& f : q_hull.facets) {
            for(size_t j = 0; j < 3; j++) {
                size_t a = f.vertices[j];
                size_t b = f.vertices[(j+1)%3];
                acc[a] += v[b];
                cnt[a]++;
            }
        }
        for(size_t j = 0; j < v.size(); j++) {
            v[j] = v[j] * 0.1 + acc[j] / cnt[j] * 0.9;
        }
        for(auto & p : v) {
            p = p * (1024 / length(p));
        }
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json("E:\\dev\\js\\qhull-js\\src\\data.js");
}
void test_3_v() {
    std::vector<vx::Vec<long double, 3>> v;
    for (size_t i = 0; i < 1024; i++) {
        double x = rand() % 2001 - 1000;
        double y = rand() % 2001 - 1000;
        if (x*x+y*y > 1000*1000) continue;
        v.emplace_back();
        v.back()[0] = x;
        v.back()[1] = y;
        v.back()[2] = (x*x+y*y)/1000;
        std::cout << v.back() << std::endl;
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json("E:\\dev\\js\\qhull-js\\src\\data.js");
}
void test_4_v() {
    std::vector<vx::Vec<long double, 4>> v;
    for (size_t i = 0; i < 8; i++) {
        double x = rand() % 2001 - 1000;
        double y = rand() % 2001 - 1000;
        double z = rand() % 2001 - 1000;
        v.emplace_back();
        v.back()[0] = x;
        v.back()[1] = y;
        v.back()[2] = z;
        v.back()[3] = (x*x+y*y+z*z)/1000;
        std::cout << v.back() << std::endl;
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json("E:\\dev\\js\\qhull-js\\src\\data.js");
}

void test_2_r() {
    std::vector<vx::Vec<long double, 2>> v;
    size_t n = 1028 * 2 * 4;
    for (size_t i = 0; i < n; i++) {
        double a = i * 3.14 * 2 / n;
        double dx = (rand()%100) / 10000.0;
        double dy = (rand()%100) / 10000.0;
        v.emplace_back(cos(a) * 1024, sin(a) * 102);
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    //test_2();
    test_3_v();
    std::cout << "DONE" << std::endl;
    //test_2_r();
}