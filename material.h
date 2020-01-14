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
#include "hitable_list.h"

static const int DIFFUSE = 1;
static const int PHONG_SPECULAR = 2;
static const int SPECULAR = 3;
static const int REFRACTION = 4;
static const int ABSORTION = 0;

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
public:
    vec4 Kd = vec4();
    vec4 Ks = vec4();
    vec4 Ksp = vec4();
    vec4 Kr = vec4();
    float ref_idx=0;
    float alpha = 0;

    bool isLight = false;
    float intensity = 0;
    vec4 color = vec4();

    material() {};
};
bool refract2(vec4 v,  vec4 n,float ref, vec4& refracted) {
    float dt= dot(normalize(v),n);
    float discrim = 1.0 - ref*ref*(1-dt*dt);
    if (discrim>0){
        refracted=ref*(normalize(v)-n*dt)-n*std::sqrt(discrim);
        return true;
    }
    return false;
}

bool refract3(hit_record rec, float cosine) {
    float r0 = (1-rec.mat->ref_idx) / (1+ rec.mat->ref_idx);
    r0=r0*r0;
    return r0 + (1-r0)*std::pow((1-cosine),5);
}


vec4 refracted_vector(vec4 rec,ray rayin, float cosine) {
    /* float r0 = (1-rec.mat->ref_idx) / (1+ rec.mat->ref_idx);
     r0=r0*r0;
     return r0 + (1-r0)*std::pow((1-cosine),5);
 */

    float c1 = dot(rec , rayin.direction)*(-1);
    float c2 = sqrt( 1 - pow(cosine,2) * (1 - pow(c1,2)));
    const vec4 refracted = (cosine * rayin.direction) + (cosine * c1 - c2) * rec;

    return refracted;
}

vec4 refract(hit_record rec, ray rayin,const hitable_list &world) {



    float t_min = 0.0;
    float t_max = std::numeric_limits<float>::max(); // Max float
    const vec4 inside_ray_origin = rec.p + .0001 * rayin.direction;
    const vec4 inside_ray_dir = refracted_vector(rec.normal,rayin, 1/rec.mat->ref_idx);
    const ray inside_ray = ray(inside_ray_origin, inside_ray_dir);
    world.hit(rayin, t_min, t_max, rec);

    // Refracted ray outside
    const vec4 out_ray_origin = rec.p ;
    const vec4 out_ray_dir = refracted_vector(rec.normal*-1,rayin, rec.mat->ref_idx);
    const ray out_ray = ray(out_ray_origin, out_ray_dir);
    return out_ray_dir;
    //const hit_record sr_out = sr.world_ptr->hit_bare_bones_obj(out_ray, std::vector<GeometryObject*>{sr.obj_ptr});

    /*  vec4 outward_normal;
      vec4 reflected = rayin.direction - 2 * rec.normal * dot(rayin.direction, rec.normal);
      float ni_over_nt;
      vec4 attenuation = vec4(1, 1, 1, 0);
      vec4 refracted;
      float reflect_prob;
      float cosine;
      if (dot(rayin.direction, rec.normal) > 0) {
          outward_normal = rec.normal*(-1);
          ni_over_nt = rec.mat->ref_idx;
          cosine = rec.mat->ref_idx*dot(rayin.direction,rec.normal)/modulus(rayin.direction);
      }
      else{
          outward_normal = rec.normal;
          ni_over_nt = 1.0/rec.mat->ref_idx;
          cosine = -dot(rayin.direction,rec.normal)/modulus(rayin.direction);
      }
      if(refract2(rayin.direction,outward_normal,ni_over_nt,refracted)){
          reflect_prob=0.6;//=refract3(rec,cosine);
      }
      else{
          reflect_prob=1.0;
          refracted=reflected;
      }
      float e=dis(generator);
      //std::cout<< e << "comparada con" << reflect_prob << std::endl;
      if(e<0.6){
          refracted=reflected;
      }
      else{
       //   std::cout << "aaaaaa" << std::endl;
          refracted=refracted;
      }
      return refracted;*/
    /*    vec4 outward_normal;
    vec4 reflected = rayin.direction - 2 * rec.normal * dot(rayin.direction, rec.normal);
    float ni_over_nt;
    vec4 attenuation = vec4(1, 1, 0, 0);
    vec4 refracted;

    if (dot(rayin.direction, rec.normal) > 0) {
        outward_normal = rec.normal*(-1);
        ni_over_nt = rec.mat->ref_idx;
    }
    else{
        outward_normal = rec.normal;
        ni_over_nt = 1.0/rec.mat->ref_idx;
    }

    if(refract2(rayin.direction,outward_normal,ni_over_nt,refracted)){
        refracted=refracted;
    }
    else{
        refracted=reflected;
    }


    return refracted;
    */
}

vec4 lobe_sampling_random_direction(hit_record rec) {
    vec4 k = normalize(rec.normal);
    vec4 i;
    vec4 ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    while (dot(k, ran) == 1) {
        ran = normalize(vec4(dis(generator), dis(generator), dis(generator), 0));
    }
    i = cross(k, ran);
    vec4 j = cross(k, i);

    mat4 T = mat4(i, j, k, rec.p);

    float theta = std::pow(acos(dis(generator)), (1 / (rec.mat->alpha + 1)));

    float phi = 2 * M_PI * dis(generator);

    vec4 nextDirectionLocal = vec4(sin(theta) * cos(phi),
                                   sin(theta) * sin(phi),
                                   cos(theta),
                                   0);

    return T * nextDirectionLocal;
}
/*
vec4 reflaction(hit_record rec) {
    vec4 k= normalize(rec.normal);
    return
}*/

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

class dielectric : public material {
public:
    dielectric(vec4 kr, float ref_idx2) {
        Kr = kr;
        ref_idx=ref_idx2;
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

vec4 phong_specular_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    vec4 reflect = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->Kd;
    vec4 Ks = record.mat->Ks;
    float alpha = record.mat->alpha;
    vec4 sp = (Ks * ((alpha + 2) / (2)));
    float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / max(Ks)) * luminance * ((Ks * ((alpha + 2) / (2 /** M_PI*/))) *
             std::pow(std::abs(dot(wi, normalize(record.normal))), alpha));

}

vec4 phong_diffuse_BRDF(hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    wo = normalize(wo);
    wi = normalize(wi);
    vec4 reflect = normalize(wo - 2 * normalize(record.normal) * dot(normalize(wo), normalize(record.normal)));
    vec4 Kd = record.mat->Kd;
    vec4 Ks = record.mat->Ks;
    float alpha = record.mat->alpha;
    vec4 sp = (Ks * ((alpha + 2) / (2)));
    float pw = std::pow(std::max(dot(wi, normalize(record.normal)), (float) 0), alpha);
    return (1 / (max(Kd))) * luminance * Kd;

}

vec4 BRDF(int event, hit_record record, vec4 wo, vec4 wi, vec4 luminance) {
    switch (event) {
        case ABSORTION:
            return vec4();
        case DIFFUSE:
            return phong_BRDF(record, wo, wi, luminance);
            return phong_diffuse_BRDF(record, wo, wi, luminance);
        case PHONG_SPECULAR:
            return phong_BRDF(record, wo, wi, luminance);
            return phong_specular_BRDF(record, wo, wi, luminance);
        case SPECULAR:
            return specular_BRDF(record, wo, wi, luminance);
        case REFRACTION:
            //return (1 / max(record.mat->Kr)) * luminance * record.mat->Kr; // TODO: Falta el seno
            return (record.mat->Kr*luminance);

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
