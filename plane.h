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

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override;
};

std::ostream &operator<<(std::ostream &os, const plane &p) {
    return os << "Center: " << p.center << "\nNormal: " << p.normal;
}

bool plane::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {
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

#endif //TRABAJO1_PLANE_H
