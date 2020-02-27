cmake_minimum_required(VERSION 3.14)
project(VoronoiX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,1000000000")

set(
        HEADER_FILES
        include/vx/core.hpp
        include/vx/math/vec.hpp
        include/vx/math/mat.hpp
        include/vx/math/qhull.hpp
        include/vx/math/geometry.hpp
        include/vx/model/graph.hpp
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