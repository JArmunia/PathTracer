//
// Created by muniia on 19/09/19.
//

#ifndef TRABAJO1_SPHERE_H
#define TRABAJO1_SPHERE_H


#include <cmath>
#include <iostream>
#include "material.h"

class sphere : public hitable {


    const vec4 center;
public:
    const float radius;
    const float radius2; // radius squared
    material *mat;

    mat4 basis = identity();

    sphere(vec4 center, float radius, material *m) : center(center), radius(radius),
                                                     radius2(radius * radius) {
        center = vec4(0,0,0,1);
        basis = basis * translation(center);
        mat = m;
    };

    bool hit(const ray &loc_ray, float t_min, float t_max, hit_record &rec) override;

};

bool sphere::hit(const ray &r, float t_min, float t_max, hit_record &rec) {
    mat4 inv = inverse(basis);
    ray loc_ray(inv * r.origin, inv * r.direction);
    vec4 oc = loc_ray.origin - center;
    float a = dot(loc_ray.direction, loc_ray.direction);
    float b = dot(oc, loc_ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - a * c;
    bool intersect = false;
    vec4 p, normal;
    if (discriminant > 0) {
        float temp = (-b + std::sqrt(b * b - a * c)) / a;
        if (t_min < temp && temp < t_max) {
            rec.t = temp;
            rec.p = basis * loc_ray.intersection(rec.t);
            rec.normal = normalize((rec.p - basis * center ));
            // Para que la normal sea positiva también desde dentro de la esfera
            rec.mat = mat;
            intersect = true;
        }

        float temp2 = (-b - std::sqrt(b * b - a * c)) / a;
        if (temp2 < temp && t_min < temp2 && temp2 < t_max) { // Podria petar aqui
            rec.t = temp2;
            rec.p = basis * loc_ray.intersection(rec.t);
            rec.normal = normalize((rec.p - basis * center));
            rec.mat = mat;

            intersect = true;
        }
        vec4 dir = normalize((inv * rec.p) - vec4(0, 0, 0, 1));
        float phi = atan2(dir.x(), dir.z());
        float theta = asin(dir.y());

        rec.v = 0.5 + (theta / M_PI);
        rec.u = (phi + M_PI) / (M_PI * 2);


    }

    return intersect;
}

#endif //TRABAJO1_SPHERE_H
