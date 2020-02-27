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
        double a = rand() % 200;
        double b = rand() % 150;
        double y = sin(a) * 512;
        double x = cos(b) * cos(a) * 512;
        double z = sin(b) * cos(a) * 512;
        v.emplace_back(x, y, z);
        std::cout << v.back() << std::endl;
    }
    vx::QHull q_hull(v);
    std::cout << q_hull.facets.size() << std::endl;
    q_hull.save_json("out.json");
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
    test_3_r();
    std::cout << "DONE" << std::endl;
    //test_2_r();

}