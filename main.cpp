#include <iostream>
#include <fstream>
#include <random>
#include "camera.h"
#include "sphere.h"
#include "plane.h"
#include "light.h"
#include "image.h"

//class material;
struct options {
    int resolution_x = 400;
    int resolution_y = 300;
    float h_fov = 120;

    int rays_per_pixel = 500;
    float shadow_bias = 10e-4;

} opt;

vec4 color(ray r, const hitable_list &world, const std::vector<point_light *> &lights, int n) {
    hit_record rec, rec1;
    bool hit_anything;
    bool hit_light;
    float t_min = 0.0;
    float t_max = std::numeric_limits<float>::max(); // Max float

    hit_record rec_min;
    vec4 rgb = vec4(0., 0., 0., 0);

    hit_anything = world.hit(r, t_min, t_max, rec); // Rayo desde la camara
    if (hit_anything) {
        ray scattered;
        int event = russian_roulette(rec.mat);

        for (auto light: lights) {
            hit_light = shadow_ray(rec, *light, world, opt.shadow_bias);
            if (hit_light) {
                vec4 light_direction = light->position - rec.p;

                //rgb = rgb + rec.mat->attenuation(rec, r.direction, light_direction,
                //                                 luminance(*light, ray(rec.p, light->position - rec.p)), 1);// *

                rgb = rgb + BRDF(event, rec, r.direction, light_direction,
                                 luminance(*light, ray(rec.p, light->position - rec.p)));
            } else {
                // rgb += 0
            }
        }
        if (event != ABSORTION) {
            if (event == DIFFUSE) {
                scattered = ray(rec.p, cosine_sampling_random_direction(rec));
            } else if (event == SPECULAR) {
                vec4 direction = r.direction - 2 * rec.normal * dot(r.direction, rec.normal);
                scattered = ray(rec.p, direction);
            } else if (event == PHONG_SPECULAR) {
                //vec4 direction = r.direction - 2 * rec.normal * dot(r.direction, rec.normal);
                scattered = ray(rec.p, cosine_sampling_random_direction(rec));
            }
//TODO: No le estoy diciendo la direccion del siguiente rayo
            //rgb = rgb + rec.mat->attenuation(rec, r.direction, scattered.direction,
            //                                 color(scattered, world, lights, n + 1), 1);
            rgb = rgb + BRDF(event, rec, r.direction, scattered.direction,
                             color(scattered, world, lights, n + 1)) ;
        }
    } else {
        //rgb += 0
    }

    return rgb;
}

//std::mt19937 generator(clock());
//std::uniform_real_distribution<double> dis(0.0, 1.0);
int main() {

    if (false) {
        camera cam = camera(vec4(0, 0, 0, 1), vec4(1, 0, 0, 0),
                            vec4(0, 1, 0, 0), vec4(0, 0, -1, 0),
                            opt.resolution_x, opt.resolution_y, opt.h_fov * M_PI / 180);

        // Merssene twister PRNG


        std::vector<point_light *> lights;


        lights.push_back(new point_light(vec4(0, 1., -4, 1),
                                         vec4(1, 1, 1, 0), 6000));


        hitable_list world;

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// SPHERES ///////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        world.hit_vector.push_back(new sphere(vec4(-3.5, -2.5, -8, 1), 1.5,
                                              new phong(vec4(0.5, 0.2, 0.2, 0),
                                                        vec4(0.1, 0.1, 0.1, 0),
                                                        10)));
        world.hit_vector.push_back(new sphere(vec4(3.5, -2.5, -8, 1), 1.5,
                                              new specular(vec4(0.8, 0.8, 0.8, 0))));
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// PLANES ////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        world.hit_vector.push_back(new plane(vec4(0, 4, 0, 1), vec4(0, -1, 0, 0),
                                             new lambertian(vec4(0.4, 0.4, 0.4, 0)))); //Superior
        world.hit_vector.push_back(new plane(vec4(0, -4, 0, 1), vec4(0, 1, 0, 0),
                                             new lambertian(vec4(0.4, 0.4, 0.4, 0)))); //Inferior
        world.hit_vector.push_back(new plane(vec4(0, 0, -10, 1), vec4(0, 0, 1, 0),
                                             new lambertian(vec4(0.4, 0.4, 0.4, 0)))); //Frontal
        world.hit_vector.push_back(new plane(vec4(5, 0, 0, 1), vec4(-1, 0, 0, 0),
                                             new lambertian(vec4(0.1, 0.4, 0.1, 0)))); //Derecho
        world.hit_vector.push_back(new plane(vec4(-5, 0, 0, 1), vec4(1, 0, 0, 0),
                                             new lambertian(vec4(0.4, 0.1, 0.1, 0))));//Izquierdo
        //world.hit_vector.push_back(new plane(vec4(0, 0, 1, 1), vec4(0, 0, -1, 0),
        //                                     new lambertian(vec4(0.8, 0.4, 0.1, 0)))); //Trasera





        vec4 rgb;
        vec4 temp;
        vec4 point;
        std::vector<std::array<float, 3>> pixels;
        float color_resolution = 0, max_rgb;
        clock_t start, end;
        start = clock();
        for (int j = 0; j < cam.resolution_y; j++) {
            for (int i = 0; i < cam.resolution_x; i++) {
                rgb = vec4(0, 0, 0, 0);
                for (int rays = 0; rays < opt.rays_per_pixel; rays++) { // Antialiasing

                    point = vec4((-cam.plane_x_size / 2) + i * (cam.plane_x_size / cam.resolution_x) +
                                 dis(generator) * (cam.plane_x_size / cam.resolution_x),
                                 (cam.plane_y_size / 2) - j * (cam.plane_y_size / cam.resolution_y) -
                                 dis(generator) * (cam.plane_y_size / cam.resolution_y),
                                 -1,
                                 1);

                    ray r = ray(cam.origin, normalize(point - cam.origin));
                    hit_record rec;

                    int o = 0;
                    int num = 0;
                    temp = color(r, world, lights, 0) / (float) opt.rays_per_pixel;
                    rgb = rgb + temp;// (float) rays_per_pixel;

                }

                max_rgb = std::max(rgb.r(), rgb.b());
                max_rgb = std::max(max_rgb, rgb.g());
                if (max_rgb > color_resolution) {
                    color_resolution = max_rgb;
                }
                pixels.push_back({rgb.r(), rgb.g(), rgb.b()});
            }
            std::cout << "\r" << 100 * j / cam.resolution_y << "%";
        }
        end = clock();
        std::cout << "\nExecution time: " << (double) (end - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

        //myfile.close();
        image hdr = image("P3", cam.resolution_x, cam.resolution_y, color_resolution, pixels);
        hdr.save("image_hdr.ppm");
        image eq = equalize(hdr, 1023);
        eq.save("image_eq.ppm");
        image gamm = equalize(gamma(hdr, 1.5), 1023);
        gamm.save("image_gamma.ppm");
        image cl = equalize(clamp(hdr, color_resolution * 0.8), 1023);
        cl.save("image_cl.ppm");
        image i = gamma(clamp(hdr, 10000), 2);
        i.save("image_fix.ppm");
        image i2 = gamma(clamp(hdr, 6000), 2);
        i2.save("image_fix2.ppm");
        image i3 = equalize(gamma(clamp(hdr, 6000), 2),1023);
        i3.save("image_fix3.ppm");
    } else {
        image i = image("image_hdr.ppm");
        i = gamma(clamp(i, 1000), 2.2);
        i.save("image_fix4.ppm");
    }
    return 0;
}
