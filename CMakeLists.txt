cmake_minimum_required(VERSION 3.1)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(Slider)

add_compile_definitions(LOCAL)
add_executable(Slider Slider.cpp)
