//
// Created by muniia on 29/10/19.
//

#ifndef TRABAJO3_CAMERA_H
#define TRABAJO3_CAMERA_H

#include "ray.h"
#include "material.h"

class camera {

    const vec4 origin = vec4(0, 0, 0, 1);
    const vec4 x = vec4(1, 0, 0, 0);
    const vec4 y = vec4(0, 1, 0, 0);
    const vec4 z = vec4(0, 0, -1, 0);
public:
    const int resolution_x;
    const int resolution_y;

    mat4 basis = identity();
    float plane_x_size;
    float plane_y_size;

/*
 * camera(vec4(0, 0, 0, 1),
 * vec4(1, 0, 0, 0),
                            vec4(0, 1, 0, 0), vec4(0, 0, -1, 0),
                            opt.width, opt.height, opt.h_fov * M_PI / 180);
 */
    camera(int resolution_x, int resolution_y, float h_fov)
            : resolution_x(resolution_x),
              resolution_y(resolution_y) {
        plane_x_size = tan(h_fov / 2);
        float aspect_ratio = (float) resolution_x / (float) resolution_y;
        plane_y_size = plane_x_size / aspect_ratio;

    };

    ray get_ray(int i, int j) {
        vec4 point((-plane_x_size / 2) + i * (plane_x_size / resolution_x) +
                   dis(generator) * (plane_x_size / resolution_x),
                   (plane_y_size / 2) - j * (plane_y_size / resolution_y) -
                   dis(generator) * (plane_y_size / resolution_y),
                   -1,
                   1);

        ray r = ray(origin, normalize(point - origin));
        vec4 origin_global = basis * origin;
        return ray(origin_global, normalize(basis * point - origin_global));

    }

};


#endif //TRABAJO3_CAMERA_H
