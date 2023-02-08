cmake_minimum_required(VERSION 3.16)
project(cdata)

set(CMAKE_CXX_STANDARD 14)

link_libraries(gsl)
include_directories(sax/include)
add_subdirectory(sax)
add_executable(cdata main.cpp)
target_link_libraries(cdata sax)


