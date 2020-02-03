//
// Created by muniia on 19/09/19.
//

#ifndef TRABAJO1_PLANE_H
#define TRABAJO1_PLANE_H


#include <iostream>
#include "hitable.h"

class plane : public hitable {

public:
    const vec4 center;
    const vec4 normal;
    material *mat;

    plane(vec4 center, vec4 normal, material *m) : center(center), normal(normalize(normal)) {
        mat = m;
    };

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) override;
};

std::ostream &operator<<(std::ostream &os, const plane &p) {
    return os << "Center: " << p.center << "\nNormal: " << p.normal;
}

bool plane::hit(const ray &r, float t_min, float t_max, hit_record &rec) {
    /*TODO comprobar que hacemos bien este hit, he negado la normal para que podamos tener
     * que apunta hacia nosotros pero no me termina de convencer */
    vec4 n = -1 * normalize(normal);
    vec4 l = normalize(r.direction);
    float denom = dot(n, l);
    vec4 p0_l0 = center - r.origin;
    if (denom > 1e-6) {
        // std::cout << "1" << std::endl;
        float temp = dot(p0_l0, n) / denom;
        if (t_min < temp && temp < t_max) {
            //  std::cout << "2" << std::endl;

            rec.t = temp;
            rec.mat = mat;
            rec.normal = normalize(normal);
            rec.p = r.intersection(rec.t);

            return (rec.t > 0);
        }
    }

    return false;
}


class finite_plane : public hitable {
    vec4 normal = vec4(0, 0, 1, 0);
    vec4 center = vec4(0, 0, 0, 1);
    float size_u;
    float size_v;
public:
    mat4 basis = identity();
    material *mat;

    finite_plane(vec4 center, float size_u, float size_v, material *mat) {
        this->size_u = size_u;
        this->size_v = size_v;
        basis = basis * translation(center);
        this->mat = mat;
    };

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) override;
};

bool finite_plane::hit(const ray &r, float t_min, float t_max, hit_record &rec) {
    /*TODO comprobar que hacemos bien este hit, he negado la normal para que podamos tener
     * que apunta hacia nosotros pero no me termina de convencer */
    mat4 inv = inverse(basis);
    ray loc_ray(inv * r.origin, inv * r.direction);
    vec4 n = -1 * normalize(normal);
    vec4 l = normalize(loc_ray.direction);
    float denom = dot(n, l);
    vec4 p0_l0 = center - loc_ray.origin;
    if (denom > 1e-6) {
        // std::cout << "1" << std::endl;
        float temp = dot(p0_l0, n) / denom;
        if (t_min < temp && temp < t_max) {
            //  std::cout << "2" << std::endl;
            vec4 intersection_point = loc_ray.intersection(temp);

            if (intersection_point.x() < -size_u / 2 || intersection_point.x() > size_u / 2 ||
                intersection_point.y() < -size_v / 2 || intersection_point.y() > size_v / 2) {
                return false;
            } else {

                rec.t = temp;
                rec.p = basis * intersection_point;
                rec.u = 0.5 + (intersection_point.x() / size_u);
                rec.v = 0.5 + (intersection_point.y() / size_v);


                rec.mat = mat;
                rec.normal = normalize(basis * normal);

                return (rec.t > 0);
            }
        }
    }

    return false;
}

#endif //TRABAJO1_PLANE_H
