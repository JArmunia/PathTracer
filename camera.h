//
// Created by muniia on 29/10/19.
//

#ifndef TRABAJO3_CAMERA_H
#define TRABAJO3_CAMERA_H

#include "ray.h"


class camera {
public:
    const vec4 origin;
    const vec4 x;
    const vec4 y;
    const vec4 z;
    const int resolution_x;
    const int resolution_y;

    float plane_x_size;
    float plane_y_size;


    camera(vec4 origin, vec4 x, vec4 y, vec4 z, int resolution_x, int resolution_y, float h_fov)
            : origin(origin), x(x), y(y), z(z),
              resolution_x(resolution_x),
              resolution_y(resolution_y) {
        plane_x_size = tan(h_fov/2);
        float aspect_ratio = (float)resolution_x / (float)resolution_y;
        plane_y_size = plane_x_size / aspect_ratio;

    };

};


#endif //TRABAJO3_CAMERA_H
