cmake_minimum_required(VERSION 3.0)
project(MyProject)

add_executable(execute algoritam.cpp)

find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(gmp REQUIRED)

target_link_libraries(execute OpenGL::GL glfw gmp)