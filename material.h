//
// Created by Javie on 01/12/2019.
//

#ifndef TRABAJO3_MATERIAL_H
#define TRABAJO3_MATERIAL_H

//#include "vec4.h"

#include <array>
#include <ctime>
//#include "texture.h"
#include "hitable.h"

static const int DIFFUSE = 1;
static const int PHONG_SPECULAR = 2;
static const int SPECULAR = 3;
static const int REFRACTION = 4;
static const int ABSORTION = 0;

std::mt19937 generator(1);
std::uniform_real_distribution<double> dis(0.0, 1.0);

vec4 cosine_sampling_random_direction(hit_record rec) {
    vec4 k = normalize(rec.normal);
    vec4 i;
    vec4 ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    while (dot(k, ran) == 1) {
        ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    }
    i = cross(k, ran);
    vec4 j = cross(k, i);

    mat4 T = mat4(i, j, k, rec.p);

    float theta = acos(sqrt(1 - dis(generator)));

    float phi = 2 * M_PI * dis(generator);

    vec4 nextDirectionLocal = vec4(sin(theta) * cos(phi),
                                   sin(theta) * sin(phi),
                                   cos(theta),
                                   0);

    return T * nextDirectionLocal;
}


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
    vec4 Kd = vec4();
    vec4 Ks = vec4();
    vec4 Ksp = vec4();
    vec4 Kr = vec4();
    float alpha = 0;


    material() {};

    //virtual vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) = 0;
//
    //virtual int scatter(const ray &ray_in, const hit_record &rec, ray &scattered) const {
    //    return 1;
    //};

};

int russian_roulette(material *mat) {
    float rr = dis(generator);
    float Kd = max(mat->Kd);
    float Ks = max(mat->Ks);
    float Ksp = max(mat->Ksp);
    float Kr = max(mat->Kr);
    if (rr > Kd + Ks + Ksp + Kr) {
        return ABSORTION;
    } else if (rr > Kd + Ks + Ksp) {
        return REFRACTION;
    } else if (rr > Kd + Ks) {
        return SPECULAR;
    } else if (rr > Kd) {
        return PHONG_SPECULAR;
    } else {
        return DIFFUSE;
    }
}

class lambertian : public material {
public:
    //vec4 Kd; // Kd

    float maxKd;

    //virtual inline int type() { return DIFFUSE; };
//
    lambertian(vec4 kd) {
        Kd = kd;
        maxKd = max(Kd);
    }
//
    //int scatter(const ray &ray_in, const hit_record &rec,
    //            ray &scattered) const override {
//
    //    float rr = dis(generator);
//
    //    // Eje de coordenadas local
//
    //    vec4 nextDirectionGlobal = cosine_sampling_random_direction(rec);
    //    vec4 randa = random_in_unit_sphere();
    //    scattered = ray(rec.p, nextDirectionGlobal);
    //    if (rr >= maxKd) {
    //        return ABSORTION;
    //    } else {
    //        return DIFFUSE;
    //    }
    //}
//
    //vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) override {
    //    return (1 / maxKd) * luminance * Kd;// M_PI *
    //    //std::max(dot(normalize(x.normal), normalize(wi)), 0.f);
    //}
};

vec4 diffuse_BRDF(hit_record x, vec4 wo, vec4 wi, vec4 luminance) {
    return (1 / max(x.mat->Kd)) * luminance * x.mat->Kd;
}

class phong : public material {
public:




    //virtual inline int type() { return PHONG; };

    phong(vec4 kd, vec4 ks, float a) {
        Kd = kd;
        Ks = ks;
        alpha = a;
    }

    //int scatter(const ray &ray_in, const hit_record &rec,
    //            ray &scattered) const override {
    //    float rr = dis(generator);
//
    //    scattered = ray(rec.p, random_in_unit_sphere());
    //    ray reflected = ray(rec.p, -1 * ray_in.direction - 2 * rec.normal
    //                                                       * dot(ray_in.direction,
    //                                                             rec.normal)); // TODO: arreglar, poner el rayo reflejado
    //    if (rr >= maxKd + maxKs) {
    //        return ABSORTION;
    //    } else if (rr >= maxKd) {
    //        vec4 direction = ray_in.direction - 2 * rec.normal * dot(ray_in.direction, rec.normal);
    //        scattered = ray(rec.p, direction);
    //        return SPECULAR;
    //    } else {
    //        scattered = ray(rec.p, random_in_unit_sphere());
    //        return DIFFUSE;
    //    }
    //}
//
    //vec4 attenuation(hit_record x, vec4 wi, vec4 wo, vec4 luminance, int event) override {
//
    //    ray reflected = ray(x.p, wo);
    //    (1 / maxKd) * luminance * diffuse_coefficient;
    //    vec4 sp = specular_coefficient * ((shininess + 2) / (2));
    //    float doot = std::abs(dot(reflected.direction, wi));
    //    float pw = pow(std::abs(dot(reflected.direction, wi)), shininess);
    //    return (1 / (maxKd + maxKs)) * luminance * ((diffuse_coefficient) +
    //                                                ((specular_coefficient * ((shininess + 2) / (2))) *
    //                                                 pow(std::abs(dot(reflected.direction, wi)), shininess)));
    //}


};

class specular : public material {
public:

    virtual inline int type() { return SPECULAR; };

    specular(vec4 ks) {
        Ksp = ks;
    }

    //int scatter(const ray &ray_in, const hit_record &rec,
    //            ray &scattered) const override {
//
    //    float rr = dis(generator);
//
    //    vec4 direction = ray_in.direction - 2 * rec.normal * dot(ray_in.direction, rec.normal);
    //    scattered = ray(rec.p, direction);
    //    //std::cout << dot(normalize(ray_in.direction),
    //    //                 normalize(rec.normal)) << std::endl;
//
    //    //if (rr >= maxKs) {
    //    //    return ABSORTION;
    //    //} else {
    //    //    return SPECULAR;
    //    //}
    //}
//
    //vec4 attenuation(hit_record x, vec4 wo, vec4 wi, vec4 luminance, int event) override {
    //    vec4 reflect = wo - 2 * normalize(x.normal) * dot(normalize(wo), normalize(x.normal));
    //    //if (wi == reflect) {
    //    //    return (1 / maxKs) * luminance *
    //    //           Ksp;/// (M_PI);//* dot(normalize(wi), normalize(x.normal)));
    //    //} else {
    //    //    return vec4(0, 0, 0, 0);
    //    //}
    //}

};

vec4 specular_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    vec4 reflect = wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal));
    if (wi == reflect) {
        return (1 / max(record.mat->Ksp)) * luminance * record.mat->Ksp; // TODO: Falta el seno
    } else {
        return vec4(0, 0, 0, 0);
    }
}

vec4 phong_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    vec4 reflect = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->Kd;
    vec4 Ks = record.mat->Ks;
    float alpha = record.mat->alpha;
    vec4 sp = (Ks * ((alpha + 2) / (2)));
    float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / (max(Kd) + max(Ks))) * luminance * (Kd +
                                                    ((Ks * ((alpha + 2) / (2))) *
                                                     std::pow(std::abs(dot(wi, normalize(record.normal))), alpha)));

}

vec4 BRDF(int event, hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    switch (event) {
        case ABSORTION:
            return vec4();
            break;
        case DIFFUSE:
            return phong_BRDF(record, wo, wi, luminance);
        case PHONG_SPECULAR:
            return phong_BRDF(record, wo, wi, luminance);
        case SPECULAR:
            return specular_BRDF(record, wo, wi, luminance);
            //case REFRACTION:
            //    return diffuse_BRDF(mat, x, wo, wi, luminance);
        default:
            std::cout << "BRDF fail" << std::endl;
            return vec4();
    }
}

//std::ostream &operator<<(std::ostream &os, const material &m) {
//    return os << "Color: R:" << m.albedo.color.r()
//              << " G: " << m.color.g()
//              << " B: " << m.color.b() << std::endl;
//}

#endif //TRABAJO3_MATERIAL_H
