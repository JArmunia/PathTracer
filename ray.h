//
// Created by muniia on 29/10/19.
//

#ifndef TRABAJO3_RAY_H
#define TRABAJO3_RAY_H


#include "algebra.h"

class ray{
public:
    vec4 origin;
    vec4 direction;

    ray(){};

    ray(const vec4 origin, const vec4 direction): origin(origin), direction(normalize(direction)){}
    vec4 intersection(float t) const {
        return origin + t * normalize(direction);
    }
};
#endif //TRABAJO3_RAY_H

