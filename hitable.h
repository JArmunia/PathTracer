//
// Created by muniia on 6/11/19.
//

#ifndef TRABAJO3_HITABLE_H
#define TRABAJO3_HITABLE_H

#include "ray.h"
class material;

struct hit_record{
    hit_record() {

    }

    float t;
    vec4 p;
    vec4 normal;
    material *mat;
};

class hitable {
public:

    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
    //virtual float hit(const vec4& center, float radius, const ray& r) const = 0;
};

#endif //TRABAJO3_HITABLE_H
