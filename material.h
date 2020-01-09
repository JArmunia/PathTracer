//
// Created by Javie on 01/12/2019.
//

#ifndef TRABAJO3_MATERIAL_H
#define TRABAJO3_MATERIAL_H

//#include "vec4.h"

#include <array>
//#include "texture.h"
#include "hitable.h"
//#include "matrix.h"


static const int LAMBERTIAN = 1;
static const int PHONG = 2;
static const int SPECULAR = 3;

static const int DIFFUSE = 1;
//static const int SPECULAR = 2;
static const int REFRACTION = 3;
static const int DEATH = 0;

std::mt19937 generator(rand()); // TODO: Poner una semilla que cambie
std::uniform_real_distribution<double> dis(0.0, 1.0);

vec4 random_in_unit_sphere() {
    static std::mt19937 generator(rand()); // TODO: Poner una semilla que cambie
    static std::uniform_real_distribution<double> dis(0.0, 2.0);
    vec4 p;
    do {
        p = vec4(dis(generator), dis(generator), dis(generator), 0) - vec4(1, 1, 1, 0);
    } while (modulus(p) >= 1.);
    return p;
}

class material {
public:

    //texture *albedo;

    //material(texture *a) : albedo(a) {}

    virtual inline int type() { return 0; };

    material() {};

    virtual vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) = 0;

    virtual int scatter(const ray &ray_in, const hit_record &rec, ray &scattered) const {
        return 1;
    };

};

class lambertian : public material {
public:
    vec4 diffuse_coefficient; // Kd
    float maxKd;

    virtual inline int type() { return LAMBERTIAN; };

    lambertian(vec4 kd) {
        diffuse_coefficient = kd;
        maxKd = max(diffuse_coefficient);
    }

    int scatter(const ray &ray_in, const hit_record &rec,
                ray &scattered) const override {

        float rr = dis(generator);

        // Eje de coordenadas local

        vec4 j = normalize(rec.normal);
        vec4 i;
        if (dot(j, vec4(1, 0, 0, 0)) != 0) {
            i = cross(j, vec4(1, 0, 0, 0));
        } else {
            i = cross(j, vec4(0, 1, 0, 0));
        }
        vec4 k = cross(i, j);

        //mat4 x = mat4(1, 2, 3, 4,
        //              5, 6, 7, 8,
        //              9, 1, 2, 3,
        //              4, 5, 6, 7);

        float theta = acos(sqrt(1 - dis(generator)));

        float phi = 2 * M_PI * dis(generator);

        vec4 nextDirectionLocal = vec4(sin(theta) * cos(phi),
                                       sin(theta) * cos(theta),
                                       cos(theta),
                                       0);

        vec4 nextDirectionGlobal;


                scattered = ray(rec.p, random_in_unit_sphere());
        if (rr >= maxKd) {
            return DEATH;
        } else {
            return DIFFUSE;
        }
    }

    vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) override {
        return (1 / maxKd) * luminance * diffuse_coefficient;// M_PI *
        //std::max(dot(normalize(x.normal), normalize(wi)), 0.f);
    }
};

class phong : public material {
public:
    vec4 diffuse_coefficient; // Kd
    vec4 specular_coefficient; //Ks
    float shininess; // alpha

    float maxKd;
    float maxKs;

    virtual inline int type() { return PHONG; };

    phong(vec4 kd, vec4 ks, float alpha) {
        diffuse_coefficient = kd;
        specular_coefficient = ks;
        shininess = alpha;
        maxKd = max(diffuse_coefficient);
        maxKs = max(specular_coefficient);
    }

    int scatter(const ray &ray_in, const hit_record &rec,
                ray &scattered) const override {
        float rr = dis(generator);

        scattered = ray(rec.p, random_in_unit_sphere());
        ray reflected = ray(rec.p, -1 * ray_in.direction - 2 * rec.normal
                                                           * dot(ray_in.direction,
                                                                 rec.normal)); // TODO: arreglar, poner el rayo reflejado
        if (rr >= maxKd + maxKs) {
            return DEATH;
        } else if (rr >= maxKd) {
            vec4 direction = ray_in.direction - 2 * rec.normal * dot(ray_in.direction, rec.normal);
            scattered = ray(rec.p, direction);
            return SPECULAR;
        } else {
            scattered = ray(rec.p, random_in_unit_sphere());
            return DIFFUSE;
        }
    }

    vec4 attenuation(hit_record x, vec4 wi, vec4 wo, vec4 luminance, int event) override {

        ray reflected = ray(x.p, wo);
        (1 / maxKd) * luminance * diffuse_coefficient;
        vec4 sp = specular_coefficient * ((shininess + 2) / (2));
        float doot = std::abs(dot(reflected.direction, wi));
        float pw = pow(std::abs(dot(reflected.direction, wi)), shininess);
        return (1 / (maxKd + maxKs)) * luminance * ((diffuse_coefficient) +
                                                    ((specular_coefficient * ((shininess + 2) / (2))) *
                                                     pow(std::abs(dot(reflected.direction, wi)), shininess)));
    }


};

class specular : public material {
public:
    vec4 specular_coefficient;
    float maxKs;

    virtual inline int type() { return SPECULAR; };

    specular(vec4 ks) {
        specular_coefficient = ks;
        maxKs = max(specular_coefficient);
    }

    int scatter(const ray &ray_in, const hit_record &rec,
                ray &scattered) const override {

        float rr = dis(generator);

        vec4 direction = ray_in.direction - 2 * rec.normal * dot(ray_in.direction, rec.normal);
        scattered = ray(rec.p, direction);
        //std::cout << dot(normalize(ray_in.direction),
        //                 normalize(rec.normal)) << std::endl;

        if (rr >= maxKs) {
            return DEATH;
        } else {
            return SPECULAR;
        }
    }

    vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) override {
        vec4 reflect = wo - 2 * normalize(x.normal) * dot(normalize(wo), normalize(x.normal));
        if (wi == reflect) {
            return (1 / maxKs) * luminance *
                   specular_coefficient;/// (M_PI);//* dot(normalize(wi), normalize(x.normal)));
        } else {
            return vec4(0, 0, 0, 0);
        }
    }

};



//std::ostream &operator<<(std::ostream &os, const material &m) {
//    return os << "Color: R:" << m.albedo.color.r()
//              << " G: " << m.color.g()
//              << " B: " << m.color.b() << std::endl;
//}

#endif //TRABAJO3_MATERIAL_H
