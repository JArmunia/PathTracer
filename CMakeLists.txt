cmake_minimum_required(VERSION 3.15)
project(PathTracer)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(PathTracer
        algebra.h
        camera.h
        hitable.h
        hitable_list.h
        image.h
        light.h
        main.cpp
        material.h
        plane.h
        ray.h
        sphere.h
        triangle.h rectangle.h)
