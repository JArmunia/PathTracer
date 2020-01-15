//
// Created by Javie on 25/11/2019.
//

#ifndef TRABAJO3_HITABLE_LIST_H
#define TRABAJO3_HITABLE_LIST_H

#include <vector>
#include "hitable.h"

class hitable_list : public hitable {
public:
    std::vector<hitable *> hit_vector;
    float intersection_bias = 0;
    hitable_list() = default;
    hitable_list(float i_bias){
        intersection_bias = i_bias;
    }
    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override;
};

bool hitable_list::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    float closest_so_far = t_max;
    // Recorremos el vector
    for (auto hittable_object : hit_vector) {
        // Comprobamos si intersecta algun objeto
        if (hittable_object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            rec = temp_rec;
            rec.p = rec.p ;
            closest_so_far = temp_rec.t;
        }
    }
    //std::cout << rec.mat << std::endl;
    return hit_anything;
}

#endif //TRABAJO3_HITABLE_LIST_H
