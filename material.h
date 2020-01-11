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

    bool isLight = false;
    float intensity = 0;
    vec4 color = vec4();

    material() {};
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

class light : public material {
public :
    light(vec4 color, float intensity) {
        this->isLight = true;
        this->color = intensity * color;
        this->intensity = intensity;
    }
};

class lambertian : public material {
public:
    lambertian(vec4 kd) {
        Kd = kd;
    }
};


class phong : public material {
public:
    phong(vec4 kd, vec4 ks, float a) {
        Kd = kd;
        Ks = ks;
        alpha = a;
    }
};

class specular : public material {
public:
    specular(vec4 ks) {
        Ksp = ks;
    }
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
