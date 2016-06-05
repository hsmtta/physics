

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(Lotkavolterra)

include_directories("C:/Program Files (x86)/GnuWin32/include")
link_directories("C:/Program Files (x86)/GnuWin32/lib")

add_executable(lotkavolterra src/main.cpp)

target_link_libraries(lotkavolterra gsl)