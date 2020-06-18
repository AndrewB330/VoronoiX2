cmake_minimum_required(VERSION 3.12)
project(VoronoiX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,1000000000")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(
        HEADER_FILES
        include/vx/core.hpp
        include/vx/geometry/vector.hpp
        include/vx/geometry/matrix.hpp
        include/vx/geometry/polyhedron.hpp
        include/vx/model/graph.hpp
        include/vx/algo/alg_convex_hull.hpp
        include/vx/algo/alg_convex_hull_2d.hpp
        include/vx/algo/alg_convex_hull_kd.hpp
        include/vx/algo/alg_delaunay_2d.hpp
        include/vx/algo/graph_delaunay.hpp
        include/vx/algo/alg_voronoi_2d.hpp
        include/vx/algo/graph_voronoi.hpp
)

set(
        SOURCE_FILES
        src/core.cpp
)

# VoronoiX Library

add_library(VoronoiX ${HEADER_FILES} ${SOURCE_FILES})

target_include_directories(
        VoronoiX
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        src
)

# VoronoiX runnable

add_executable(main main.cpp)

target_link_libraries(main VoronoiX)