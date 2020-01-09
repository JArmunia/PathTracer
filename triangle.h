//
// Created by muniia on 19/09/19.
//

#ifndef TRABAJO1_TRIANGLE_H
#define TRABAJO1_TRIANGLE_H


#include <iostream>
#include "hitable.h"

class triangle : public hitable {

public:
    const vec4 p0;
    const vec4 p1;
    const vec4 p2;
    material *mat;

    triangle(vec4 p0, vec4 p1, vec4 p2, material *m) : p0(p0),  p1(p1),  p2(p2){
        mat = m;
    };

    float intersection(ray r);

    bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override;
};

// TODO preguntar si esto se hace asi

/*std::ostream &operator<<(std::ostream &os, const triangle &s) {
    return os << "Center: " << s.center << "\nRadius: " << s.radius;
}*/

float triangle::intersection(ray r) {
  /*  vec4 oc = r.origin - center;
    float a = dot(r.direction, r.direction);
    float b = 2.0 * dot(oc, r.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return -1;
    } else {
        std::cout << (-b - sqrt(discriminant)) / (2.0 * a) << std::endl;
        return (-b - sqrt(discriminant)) / (2.0 * a);
    }*/
  return 0;
}

bool triangle::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {


   // (
            // const Vec3f &orig, const Vec3f &dir,
  //  const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
  //  float &t)
    //{
        // compute plane's normal
        vec4 p0p1 = p1 - p0;
        vec4 p0p2 = p2 - p0;
        // no need to normalize
        vec4 N = cross(p0p1,p0p2); // N
        rec.normal=N;
        rec.mat=mat;

        float area2 = modulus(N);

        // Step 1: finding P

        // check if ray and plane are parallel ?
        float NdotRayDirection = dot(N,r.direction);
        if (fabs(NdotRayDirection) < 10e-3) // almost 0
            return false; // they are parallel so they don't intersect !

        // compute d parameter using equation 2
        float d = dot(N,p0);

        // compute t (equation 3)
        double t = (dot(N,r.origin) + d) / NdotRayDirection;
        // check if the triangle is in behind the ray
        if (t < 0) return false; // the triangle is behind
        if (t_min > t && t > t_max) {
            return false;
        }
        // compute the intersection point using equation 1
        vec4 P = r.origin + t * r.direction;
        rec.p=P;
        rec.t=t;
        // Step 2: inside-outside test
        vec4 C; // vector perpendicular to triangle's plane

        // edge 0
        vec4 edge0 = p1 - p0;
        vec4 pp0 = P - p0;
        C = cross(edge0,pp0);
        if (dot(N,C) < 0) return false; // P is on the right side

        // edge 1
        vec4 edge1 = p2 - p1;
        vec4 pp1 = P - p1;
        C = cross(edge1,pp1);
        if (dot(N,C) < 0)  return false; // P is on the right side

        // edge 2
        vec4 edge2 = p0 - p2;
        vec4 pp2 = P - p2;
        C = cross(edge2,pp2);
        if (dot(N,C) < 0) return false; // P is on the right side;
        return true;
/*    Vector    u, v, n;             // triangle vectors
    Vector    dir, w0, w;          // ray vectors
    float     r, a, b;             // params to calc ray-plane intersect

    // get triangle edge vectors and plane normal
    u = p1 - p0;
    v = p2 - p0;
    n = u * v;             // cross product
    if (n == [0,0,0])            // triangle is degenerate
        return false;                 // do not deal with this case

    dir = p1 - p0;             // ray direction vector
    w0 = p0 - p0;
    a = -dot(n,w0);
    b = dot(n,dir);
    if (fabs(b) < 0.00000001) {     // ray is parallel to triangle plane
        return false;             // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
        return false;                  // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    *I = R.P0 + r * dir;           // intersect point of ray and plane

    // is I inside T?
    float    uu, uv, vv, wu, wv, D;
    uu = dot(u,u);
    uv = dot(u,v);
    vv = dot(v,v);
    w = *I - T.V0;
    wu = dot(w,u);
    wv = dot(w,v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    float s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
        return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return false;

    return true;                      // I is in T
    return intersect;
    */
}

#endif //TRABAJO1_SPHERE_H
