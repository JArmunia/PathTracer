//
// Created by muniia on 19/09/19.
//

#ifndef TRABAJO1_SPHERE_H
#define TRABAJO1_SPHERE_H


#include <iostream>
#include "material.h"

class sphere : public hitable {

public:
    const vec4 center;
    const float radius;
    const float radius2; // radius squared
    material *mat;

    sphere(vec4 center, float radius, material *m) : center(center), radius(radius),
                                                    radius2(radius * radius) {
        mat = m;
    };

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override;
};

std::ostream &operator<<(std::ostream &os, const sphere &s) {
    return os << "Center: " << s.center << "\nRadius: " << s.radius;
}

bool sphere::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {
    vec4 oc = r.origin - center;
    float a = dot(r.direction, r.direction);
    float b = dot(oc, r.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - a * c;
    bool intersect = false;
    vec4 p, normal;
    if (discriminant > 0) {
        float temp = (-b + sqrt(b * b - a * c)) / a;
        if (t_min < temp && temp < t_max) {
            rec.t = temp;
            rec.p = r.intersection(rec.t);
            rec.normal = (rec.p - center) / radius;
            // Para que la normal sea positiva también desde dentro de la esfera
            rec.normal = -sgn(dot(rec.normal, r.direction)) * rec.normal;
            rec.mat = mat;
            //std::cout << mat << std::endl;
            intersect =  true;
        }

        float temp2 = (-b - sqrt(b * b - a * c)) / a;
        if (temp2 < temp && t_min < temp2 && temp2 < t_max) { // Podria petar aqui
            rec.t = temp2;
            rec.p = r.intersection(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.normal = -sgn(dot(rec.normal, r.direction)) * rec.normal;
            rec.mat = mat;

            intersect =  true;
        }

    }

    return intersect;
}

#endif //TRABAJO1_SPHERE_H
