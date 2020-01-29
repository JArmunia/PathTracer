//
// Created by Javie on 29/01/2020.
//

#include "hitable.h"

#ifndef PATHTRACER_RECTANGLE_H
#define PATHTRACER_RECTANGLE_H

#endif //PATHTRACER_RECTANGLE_H


class xy_rect : public hitable {

public:
    float x0, x1, y0, y1, k;
    material *mat;

    xy_rect() = default;

    xy_rect(float x0, float x1, float y0, float y1, float k, material *mat) :
            x0(x0), x1(x1), y0(y0), y1(y1), k(k), mat(mat) {};


    bool hit(const ray &r, float t_min, float t_max, hit_record &rec)  override;
};

bool xy_rect::hit(const ray &r, float t_min, float t_max, hit_record &rec)  {
    float t = (k - r.origin.z()) / r.direction.z();
    if (t < t_min || t > t_max) {
        return false;
    }
    float x = r.origin.x() + t * r.direction.x();
    float y = r.origin.y() + t * r.direction.y();
    if (x < x0 || x > x1 || y < y0 || y > y1) {
        return false;
    }
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (y - y0) / (y1 - y0);
    rec.t = t;
    rec.mat = mat;
    rec.p = r.intersection(t);
    rec.normal = vec4(0, 0, 1, 0);
    return true;
}


class xz_rect : public hitable {

public:
    float x0, x1, z0, z1, k;
    material *mat;

    xz_rect() = default;

    xz_rect(float x0, float x1, float y0, float y1, float k, material *mat) :
            x0(x0), x1(x1), z0(y0), z1(y1), k(k), mat(mat) {};


    bool hit(const ray &r, float t_min, float t_max, hit_record &rec)  override;
};

bool xz_rect::hit(const ray &r, float t_min, float t_max, hit_record &rec)  {
    float t = (k - r.origin.y()) / r.direction.y();
    if (t < t_min || t > t_max) {
        return false;
    }
    float x = r.origin.x() + t * r.direction.x();
    float z = r.origin.z() + t * r.direction.z();
    if (x < x0 || x > x1 || z < z0 || z > z1) {
        return false;
    }
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    rec.mat = mat;
    rec.p = r.intersection(t);
    rec.normal = vec4(0, 1, 0, 0);
    return true;
}


class yz_rect : public hitable {

public:
    float z0, z1, y0, y1, k;
    material *mat;

    yz_rect() = default;

    yz_rect(float y0, float y1, float z0, float z1, float k, material *mat) :
            y0(y0), y1(y1), z0(z0), z1(z1), k(k), mat(mat) {};


    bool hit(const ray &r, float t_min, float t_max, hit_record &rec)  override;
};

bool yz_rect::hit(const ray &r, float t_min, float t_max, hit_record &rec)  {
    float t = (k - r.origin.x()) / r.direction.x();
    if (t < t_min || t > t_max) {
        return false;
    }
    float y = r.origin.y() + t * r.direction.y();
    float z = r.origin.z() + t * r.direction.z();
    if (y < y0 || y > y1 || z < z0 || z > z1) {
        return false;
    }
    rec.u = (y - y0) / (y1 - y0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    rec.mat = mat;
    rec.p = r.intersection(t);
    rec.normal = vec4(1, 0, 0, 0);
    return true;
}

class box : public hitable {
public:
    xy_rect front, back;
    yz_rect right, left;
    xz_rect top, bottom;
    material *mat;
    hitable_list planes;

    box(vec4 xyz0, vec4 xyz1, material *mat) {
        front = xy_rect(xyz0.x(), xyz1.x(), xyz0.y(), xyz1.y(), xyz0.z(), mat);
        back = xy_rect(xyz0.x(), xyz1.x(), xyz0.y(), xyz1.y(), xyz1.z(), mat);
        right = yz_rect(xyz0.y(), xyz1.y(), xyz0.z(), xyz1.z(), xyz0.x(), mat);
        left = yz_rect(xyz0.y(), xyz1.y(), xyz0.z(), xyz1.z(), xyz1.x(), mat);
        bottom = xz_rect(xyz0.x(), xyz1.x(), xyz0.z(), xyz1.z(), xyz0.y(), mat);
        top = xz_rect(xyz0.x(), xyz1.x(), xyz0.z(), xyz1.z(), xyz1.y(), mat);

        planes.hit_vector.push_back(&front);
        planes.hit_vector.push_back(&back);
        planes.hit_vector.push_back(&right);
        planes.hit_vector.push_back(&left);
        planes.hit_vector.push_back(&bottom);
        planes.hit_vector.push_back(&top);

    }

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec)  override;
};

bool box::hit(const ray &r, float t_min, float t_max, hit_record &rec)  {
    return planes.hit(r, t_min, t_max, rec);
}
