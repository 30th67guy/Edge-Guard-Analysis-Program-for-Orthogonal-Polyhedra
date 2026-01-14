cmake_minimum_required(VERSION 3.20)
project(PolyhedronViewer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL CONFIG REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)

add_definitions(-DCGAL_DO_NOT_USE_GMP)
add_definitions(-DCGAL_DO_NOT_USE_MPFR)

add_executable(PolyhedronViewer src/main.cpp)
target_link_libraries(PolyhedronViewer PRIVATE CGAL::CGAL OpenGL::GL glfw)
