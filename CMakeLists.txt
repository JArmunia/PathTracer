cmake_minimum_required(VERSION 3.15)
project(Trabajo3)

set(CMAKE_CXX_STANDARD 17)
include_directories(.)

add_executable(Trabajo3
        .gitignore
        camera.h
        hitable.h
        hitable_list.h
        image.h
        light.h
        main.cpp
        material.h
        matrix.h
        plane.h
        ray.h
        sphere.h
        texture.h
        triangle.h
        algebra.h)

set_property(TARGET Trabajo3 PROPERTY CXX_STANDARD 17)
