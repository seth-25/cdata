cmake_minimum_required(VERSION 3.16)
project(cdata)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -mavx2")
link_libraries(gsl)
include_directories(sax/include)
add_subdirectory(sax)
add_executable(cdata main.cpp)
add_executable(cdata_query  main_query.cpp)
target_link_libraries(cdata sax)
target_link_libraries(cdata_query sax)


