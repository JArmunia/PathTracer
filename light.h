//
// Created by Javie on 27/11/2019.
//

#ifndef TRABAJO3_LIGHT_H
#define TRABAJO3_LIGHT_H

#include "ray.h"
#include "hitable_list.h"

class point_light {
public:
    const vec4 position;
    const vec4 color;
    const float intensity;

    point_light();

    point_light(vec4 p, vec4 c, float i) : position(p), color(c), intensity(i) {};
};

vec4 luminance(point_light light, vec4 point) {
    return light.color * light.intensity / pow(modulus(light.position - point), 2);
}

/**
 * Throws a shadow ray, we prevent shadow acne by moving the origin of the shadow ray
 * in the direction of the normal of the hit the amount specified by shadow_bias
 * @param hit_rec
 * @param light
 * @param world
 * @param shadow_bias
 * @return True if the ray doesn't hit anything between the hit and the light
 */
bool shadow_ray(hit_record hit_rec, point_light light, hitable_list world, float shadow_bias) {
    vec4 shadow_ray_origin = hit_rec.p + shadow_bias * hit_rec.normal;
    vec4 shadow_ray_direction = light.position - shadow_ray_origin;
    ray shadow_ray = ray(shadow_ray_origin, shadow_ray_direction);
    hit_record shadow_ray_rec;

    return !world.hit(shadow_ray, 0, modulus(shadow_ray_direction), shadow_ray_rec);
}


#endif //TRABAJO3_LIGHT_H
