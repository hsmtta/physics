

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(Lotkavolterra)

# include_directories("C:/Program Files (x86)/GnuWin32/include")
include_directories("/usr/local/include")

# link_directories("C:/Program Files (x86)/GnuWin32/lib")
link_directories("/usr/local/lib")

add_executable(lotkavolterra src/main.cpp)

target_link_libraries(lotkavolterra gsl)