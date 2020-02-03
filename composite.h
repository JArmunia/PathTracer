//
// Created by Javie on 03/02/2020.
//

#include "hitable.h"
#include "hitable_list.h"
#include "plane.h"

#ifndef PATHTRACER_COMPOSITE_H
#define PATHTRACER_COMPOSITE_H

#endif //PATHTRACER_COMPOSITE_H

class box : public hitable {
    hitable_list list;
public:
    mat4 basis = identity();
    finite_plane top;
    finite_plane bottom;
    finite_plane front;
    finite_plane back;
    finite_plane right;
    finite_plane left;

    box(vec4 center, float heigth, float width, float depth, material *mat) :
            top(finite_plane(vec4(0, heigth / 2, 0, 1), width, depth, mat)),
            bottom(finite_plane(vec4(0, -heigth / 2, 0, 1), width, depth, mat)),
            front(finite_plane(vec4(0, 0, 0, -depth/ 2), width, heigth, mat)),
            back(finite_plane(vec4(0, 0, depth / 2, 1), width, heigth, mat)),
            right(finite_plane(vec4(width / 2, 0, 0, 1), depth, heigth, mat)),
            left(finite_plane(vec4(-width / 2, 0, 0, 1), depth, heigth, mat)) {

        top.basis = top.basis * rotation_x(-M_PI / 2);
        list.hit_vector.push_back(&top);

        bottom.basis = bottom.basis * rotation_x(M_PI / 2);
        list.hit_vector.push_back(&bottom);


        front.basis = front.basis * rotation_y(M_PI);
        list.hit_vector.push_back(&front);

        list.hit_vector.push_back(&back);


        right.basis = right.basis * rotation_y(M_PI / 2);
        list.hit_vector.push_back(&right);

        left.basis = left.basis * rotation_y(-M_PI / 2);
        list.hit_vector.push_back(&left);

        basis = basis * translation(center);

    }

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) {
        mat4 inv = inverse(basis);
        ray r_loc = ray(inv * r.origin, inv * r.direction);
        bool hitted = list.hit(r_loc, t_min, t_max, rec);
        rec.p = basis * rec.p;
        rec.normal = basis * rec.normal;
        return hitted;
    }

    void setTopMaterial(material *m) {
        box::top.mat = m;
    }

    void setBottomMaterial(material *m) {
        box::bottom.mat = m;
    }

    void setFrontMaterial(material *m) {
        box::front.mat = m;
    }

    void setBackMaterial(material *m) {
        box::back.mat = m;
    }

    void setRightMaterial(material *m) {
        box::right.mat = m;
    }

    void setLeftMaterial(material *m) {
        box::left.mat = m;
    }


};