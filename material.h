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
#include "image.h"

static const int ABSORTION = 0;
static const int DIFFUSE = 1;
static const int PHONG_SPECULAR = 2;
static const int SPECULAR = 3;
static const int REFRACTION = 4;

class material;

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
protected:
    vec4 color = vec4();
public:
    vec4 Kd = vec4();
    vec4 Ks = vec4();
    vec4 Ksp = vec4();
    vec4 Kr = vec4();
    float alpha = 0;
    float refraction_index = 0;


    bool isLight = false;
    float intensity = 0;


    material() = default;;

    virtual vec4 albedo(hit_record rec) = 0;

    virtual vec4 getLight(hit_record rec) {
        return vec4();
    };
};


vec4 lobe_sampling_random_direction(hit_record rec, vec4 reflected) {
    vec4 wr = normalize(reflected);
    vec4 k = wr;
    vec4 i;
    vec4 ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    while (dot(k, ran) == 1) {
        ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    }
    i = cross(k, ran);
    vec4 j = cross(k, i);

    mat4 T = mat4(i, j, k, rec.p);
    vec4 nextDirectionLocal;
    do {
        float theta = std::pow(acos(dis(generator)), (1 / (rec.mat->alpha + 1)));

        float phi = 2 * M_PI * dis(generator);

        nextDirectionLocal = vec4(sin(theta) * cos(phi),
                                  sin(theta) * sin(phi),
                                  cos(theta),
                                  0);
    } while (dot(inverse(T) * normalize(rec.normal), nextDirectionLocal) < 0);
    return T * nextDirectionLocal;


}

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
        this->color = color;
        this->intensity = intensity;
    }

    vec4 albedo(hit_record rec) override {
        return vec4();
    }

    vec4 getLight(hit_record rec) {
        return intensity * color;
    }
};

//class lambertian : public material {
//public:
//    lambertian(vec4 kd) {
//        Kd = kd;
//    }
//};


class phong : public material {
public:
    phong(vec4 kd, vec4 ks, float a) {
        Kd = kd;
        Ks = ks;
        alpha = a;
    }

    vec4 albedo(hit_record rec) override {
        return Kd;
    }
};

class specular : public material {
public:
    specular(vec4 ks) {
        Ksp = ks;
    }

    vec4 albedo(hit_record rec) override {
        return vec4();
    }
};

class glass : public material {
public:
    glass(vec4 kr, float rf) {
        Kr = kr;
        refraction_index = rf;
    }

    vec4 albedo(hit_record rec) override {
        return vec4();
    }
};

vec4 reflected(vec4 direction, vec4 normal) {
    direction = normalize(direction);
    normal = normalize(normal);
    return direction - 2 * normal * dot(direction, normal);
}


vec4 refract(const vec4 &I, const vec4 &N, const float &ior) {
    float cosi = dot(normalize(I), normalize(N));
    float etai = 1, etat = ior;
    vec4 n = N;
    if (cosi < 0) { cosi = -cosi; }
    else {
        std::swap(etai, etat);
        n = -1 * N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    //return k < 0 ? vec4() : eta * I + (eta * cosi - std::sqrt(k)) * n;
    return eta * I + (eta * cosi - std::sqrt(k)) * n;
}

vec4 refracted(const vec4 wo, const hit_record rec) {
    vec4 woG = normalize(wo);
    vec4 normalG = normalize(rec.normal);
    float ior = rec.mat->refraction_index;
    if (dot(woG, normalG) < 0) {
        normalG = -1 * normalG;
        ior = 1 / ior;
    }
    float sin_theta1 = std::sqrt(1 - std::pow(dot(woG, normalG), 2));

    float theta2 = std::asin(ior * sin_theta1);
    vec4 k = normalG;
    vec4 j = cross(woG, k);
    vec4 i = cross(j, k);

    mat4 T = mat4(i, j, k, rec.p);


    vec4 refracted = vec4(sin(theta2) * cos(M_PI),
                          sin(theta2) * sin(M_PI),
                          cos(theta2),
                          0);

    return T * refracted;
}


vec4 specular_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    //vec4 reflect = wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal));
    vec4 reflect = reflected(normalize(wo), normalize(record.normal));
    if (wi == reflect) {
        return (1 / max(record.mat->Ksp)) * luminance * record.mat->Ksp;
        //*std::max(1 - std::pow(std::abs(dot(normalize(wo), normalize(record.normal))), 2),0.);
        //modulus(cross(normalize(wo), normalize(record.normal))); // TODO: Falta el seno
    } else {
        return vec4(0, 0, 0, 0);
    }
}

vec4 refraction_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    //vec4 reflect = wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal));
    vec4 refracted = refract(wo, record.normal, record.mat->refraction_index);
    //if (refracted == vec4())
    //    refracted = reflected(wo, record.normal);
    if (normalize(wi) == normalize(refracted)) {
        //std::cout <<"Init: " << refracted << " | " << wi << "end" << std::endl;
        return (1 / max(record.mat->Kr)) * luminance * record.mat->Kr;
    } else {
        //std::cout <<"Init: " << refracted << " | " << wi << "end" << std::endl;
        return vec4(0., 0, 0, 0);
    }

}

vec4 phong_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    vec4 reflect = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->albedo(record);
    vec4 Ks = record.mat->Ks;
    //std::cout << Kd << std::endl;
    float alpha = record.mat->alpha;
    vec4 sp = (Ks * ((alpha + 2) / (2)));
    float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / (max(Kd) + max(Ks))) * luminance * (Kd +
                                                    ((Ks * ((alpha + 2) / (2))) *
                                                     std::pow(std::abs(dot(wi, reflect)), alpha)));

}

vec4 phong_specular_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    vec4 wr = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->Kd;
    vec4 Ks = record.mat->Ks;
    float alpha = record.mat->alpha;
    vec4 sp = (Ks * ((alpha + 2) / (2)));
    float cos_wi = dot(wi, normalize(record.normal));
    float sin_wi = std::sqrt(1 - std::pow(cos_wi, 2));
    float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / max(Ks)) * luminance * Ks * (alpha + 2) * cos_wi * sin_wi /
           ((alpha + 1) * std::sqrt(1 - std::pow(dot(wi, wo), 2)));
    //std::abs(dot(wi,wr)) * std::sqrt());

}

vec4 phong_diffuse_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    //vec4 reflect = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->albedo(record);
    //vec4 Ks = record.mat->Ks;
    //float alpha = record.mat->alpha;
    //vec4 sp = (Ks * ((alpha + 2) / (2)));
    //float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / (max(Kd))) * luminance * Kd;

}

class texture : public material {
public:
    image img;
    int nx, ny;
    float rep_x = 1, rep_y = 1;


    texture(image i, float repetitions_x, float repetitions_y) : img(i) {
        Kd = vec4(0.9, 0.9, 0.9, 0);
        nx = i.resolution[0];
        ny = i.resolution[1];
        rep_x = repetitions_x;
        rep_y = repetitions_y;
    }

    texture(image i, float repetitions_x, float repetitions_y, float intensity) : img(i) {
        Kd = vec4(0.9, 0.9, 0.9, 0);
        nx = i.resolution[0];
        ny = i.resolution[1];
        rep_x = repetitions_x;
        rep_y = repetitions_y;
        this->isLight = true;
        this->intensity = intensity;
    }

    texture(image i) : img(i) {
        Kd = vec4(0.9, 0.9, 0.9, 0);
        nx = i.resolution[0];
        ny = i.resolution[1];
    }

    texture(image i, float intensity) : img(i) {
        Kd = vec4(0.9, 0.9, 0.9, 0);
        nx = i.resolution[0];
        ny = i.resolution[1];
        this->isLight = true;
        this->intensity = intensity;
    }


    vec4 getLight(hit_record rec) override {
        return albedo(rec);
    }

    vec4 albedo(hit_record rec) override {
        //int i1 = rec.u * (float) nx;
        //int j1 = (1 - rec.v) * (float) ny - 0.001;

        int n = rec.u / (1 / rep_x);
        float resto_u = rec.u - n * (1 / rep_x);
        n = rec.v / (1 / rep_y);
        float resto_v = rec.v - n * (1 / rep_y);

        float i_real = resto_u / (1 / rep_x) * (float) nx;
        int i1 = std::floor(i_real);
        int i2 = i1 + 1;
        float variacion_i = i_real - (float) i1;

        float j_real = 1 - resto_v / (1 / rep_y) * (float) ny - 0.001;
        int j1 = std::floor(j_real);
        int j2 = j1 + 1;
        float variacion_j = j_real - (float) j1;


        float cr = img.color_resolution;


        if (!(i1 < 0) && !(j1 < 0) && !(i2 >= nx) && !(j2 >= ny)) {
            float r = img.pixels[i1 + nx * j1][0] / img.color_resolution;
            float g = img.pixels[i1 + nx * j1][1] / img.color_resolution;
            float b = img.pixels[i1 + nx * j1][2] / img.color_resolution;

            std::array<float, 3> p1 = img.pixels[i1 + nx * j1];
            std::array<float, 3> p2 = img.pixels[i2 + nx * j1];
            std::array<float, 3> p3 = img.pixels[i1 + nx * j2];
            std::array<float, 3> p4 = img.pixels[i2 + nx * j2];
            vec4 p11(p1[0] / cr, p1[1] / cr, p1[2] / cr, 0);
            vec4 p12(p3[0] / cr, p3[1] / cr, p3[2] / cr, 0);
            vec4 p21(p2[0] / cr, p2[1] / cr, p2[2] / cr, 0);
            vec4 p22(p4[0] / cr, p4[1] / cr, p4[2] / cr, 0);

            vec4 color_total = p11 * (1 - variacion_i) * (1 - variacion_j) +
                               p12 * (1 - variacion_i) * (variacion_j) +
                               p21 * (variacion_i) * (1 - variacion_j) +
                               p22 * (variacion_i) * variacion_j;



            //if (r <= 0.001) r = 0.001;
            //if (g <= 0.001) g = 0.001;
            //if (b <= 0.001) b = 0.001;
            if (color_total.r() <= 0.001) color_total = color_total + vec4(0.01, 0, 0, 0);
            if (color_total.g() <= 0.001) color_total = color_total + vec4(0, 0.01, 0, 0);
            if (color_total.b() <= 0.001) color_total = color_total + vec4(0, 0, 0.01, 0);

            if (rec.mat->isLight) {
                //return vec4(r, g, b, 0) * intensity;
                return vec4(r, g, b, 0) * intensity;
            } else {
                //return vec4(r, g, b, 0) * max(Kd);
                return vec4(r, g, b, 0) * max(Kd);
            }
        } else {
            std::cout << "!i1 < 0: " << !(i1 < 0)
                      << " !j1 < 0: " << !(j1 < 0)
                      << " !i2 >= 0: " << !(i2 >= nx)
                      << " !j2 >= 0: " << !(j2 >= ny)
                      << std::endl;
            return vec4(0.91, 0.01, 0.01, 0.00);
        }
    }
};

vec4 BRDF(int event, const hit_record record, const vec4 wo, const vec4 wi, const vec4 luminance) {
    switch (event) {
        case ABSORTION:
            return vec4();
        case DIFFUSE:
            return phong_BRDF(record, wo, wi, luminance);
        case PHONG_SPECULAR:
            return phong_BRDF(record, wo, wi, luminance);
        case SPECULAR:
            return specular_BRDF(record, wo, wi, luminance);
        case REFRACTION:
            return refraction_BRDF(record, wo, wi, luminance);// * vec4(1,1,1,0) * 1/max(record.mat->Kr);
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
