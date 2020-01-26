#include <iostream>
#include <random>
#include <thread>
#include <atomic>
#include <future>
#include "camera.h"
#include "sphere.h"
#include "plane.h"
#include "light.h"
#include "image.h"


struct options {
    int width = 400;
    int height = 300;
    float h_fov = 120;

    int rays_per_pixel = 2000;
    float shadow_bias = 10e-4;
    float intersection_bias = 10e-4;
    int cores = 8;//std::thread::hardware_concurrency();


} opt;

vec4 color(ray r, const hitable_list &world, const std::vector<point_light *> &lights, int n) {
    hit_record rec, rec1;
    bool hit_anything;
    bool hit_light;
    float t_min = 0.0;
    float t_max = std::numeric_limits<float>::max(); // Max float

    hit_record rec_min;
    vec4 rgb = vec4();

    hit_anything = world.hit(r, t_min, t_max, rec); // Rayo desde la camara
    if (hit_anything) {
        ray scattered;
        int event = russian_roulette(rec.mat);

        for (auto light: lights) {
            hit_light = shadow_ray(rec, *light, world, opt.shadow_bias);
            if (hit_light) {
                vec4 light_direction = light->position - rec.p;

                //rgb = rgb + rec.mat->attenuation(rec, r.next_direction, light_direction,
                //                                 luminance(*light, ray(rec.p, light->position - rec.p)), 1);// *

                rgb = rgb + BRDF(event, rec, r.direction, light_direction,
                                 luminance(*light, ray(rec.p, light->position - rec.p)));
            } else {
                // rgb += 0
            }
        }
        if (rec.mat->isLight)
            return rec.mat->color;
        vec4 next_direction;
        vec4 point;
        if (event != ABSORTION) {
            switch (event) {
                case DIFFUSE: {
                    next_direction = cosine_sampling_random_direction(rec);
                    rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, next_direction);
                    break;
                }

                case PHONG_SPECULAR: {
                    //scattered = ray(rec.p, cosine_sampling_random_direction(rec));
                    //next_direction = cosine_sampling_random_direction(rec);
                    next_direction = lobe_sampling_random_direction(rec, reflected(r.direction, rec.normal));
                    rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, lobe_sampling_random_direction(rec));
                    break;
                }

                case SPECULAR: {
                    //vec4 next_direction = r.next_direction - 2 * rec.normal * dot(r.next_direction, rec.normal);
                    //scattered = ray(rec.p, reflected(r.next_direction, rec.normal));
                    next_direction = reflected(r.direction, rec.normal);
                    rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    break;
                }

                case REFRACTION: {
                    next_direction = refract(r.direction, rec.normal, rec.mat->refraction_index);
                    rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, next_direction);
                    break;
                }
            }
            //rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
            scattered = ray(rec.p, next_direction);
            rgb = rgb + BRDF(event, rec, r.direction, scattered.direction,
                             color(scattered, world, lights, n + 1));
        }

    }

    return rgb;
}

vec4 trace(camera cam, int i, int j, hitable_list scene, std::vector<point_light *> lights) {
    vec4 rgb;
    vec4 point;
    for (int rays = 0; rays < opt.rays_per_pixel; rays++) { // Antialiasing

        point = vec4((-cam.plane_x_size / 2) + i * (cam.plane_x_size / cam.resolution_x) +
                     dis(generator) * (cam.plane_x_size / cam.resolution_x),
                     (cam.plane_y_size / 2) - j * (cam.plane_y_size / cam.resolution_y) -
                     dis(generator) * (cam.plane_y_size / cam.resolution_y),
                     -1,
                     1);

        ray r = ray(cam.origin, normalize(point - cam.origin));
        hit_record rec;

        rgb = rgb + color(r, scene, lights, 0);

    }
    return rgb / opt.rays_per_pixel;
}

//std::mt19937 generator(clock());
//std::uniform_real_distribution<double> dis(0.0, 1.0);
int main() {
    if (true) {
        camera cam = camera(vec4(0, 0, 0, 1), vec4(1, 0, 0, 0),
                            vec4(0, 1, 0, 0), vec4(0, 0, -1, 0),
                            opt.width, opt.height, opt.h_fov * M_PI / 180);

        // Merssene twister PRNG


        std::vector<point_light *> lights;


        //lights.push_back(new point_light(vec4(0, 1., -4, 1),
        //                                 vec4(1, 1, 1, 0), 6000));


        hitable_list scene(opt.intersection_bias);

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// SPHERES ///////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        scene.hit_vector.push_back(new sphere(vec4(0, 1.5, -8, 1), 1.,
                                              new phong(vec4(0.5, 0.2, 0.2, 0),
                                                        vec4(0.1, 0.1, 0.1, 0),
                                                        10)));
        scene.hit_vector.push_back(new sphere(vec4(0, -1.5, -8, 1), 1.,
                                              new phong(vec4(0.2, 0.2, 0.5, 0),
                                                        vec4(0.1, 0.1, 0.1, 0),
                                                        10)));
        scene.hit_vector.push_back(new sphere(vec4(0, 0, -6, 1), 1.,
                                              new glass(vec4(0.9, 0.9, 0.9, 0), 1.8)));

        //scene.hit_vector.push_back(new sphere(vec4(1.3, -.5, -8, 1), 1.,
        //                                      new specular(vec4(0.9, 0.9, 0.9, 0))));
//
        //scene.hit_vector.push_back(new sphere(vec4(2.5, -2.5, -8, 1), 1.5,
        //                                      new specular(vec4(0.9, 0.9, 0.9, 0))));
//
        scene.hit_vector.push_back(new sphere(vec4(3, -2.5, -8, 1), 1.5,
                                              new specular(vec4(0.9, 0.9, 0.9, 0))));
//
        scene.hit_vector.push_back(new sphere(vec4(-3, -2.5, -8, 1), 1.5,
                                              new phong(vec4(0.5, 0.1, 0.1, 0),
                                                        vec4(0.2, 0.2, 0.2, 0),
                                                        15)));

        //scene.hit_vector.push_back(new sphere(vec4(0, 5, -6, 1), 1.5,
        //                                      new light(vec4(1, 1, 1, 0), 10000)));
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// PLANES ////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        scene.hit_vector.push_back(new plane(vec4(0, 4, 0, 1), vec4(0, -1, 0, 0),
                                             new light(vec4(1, 1, 1, 0), 10000))); //Superior

        scene.hit_vector.push_back(new plane(vec4(0, -4, 0, 1), vec4(0, 1, 0, 0),
                                             new phong(vec4(0.7, 0.7, 0.7, 0),
                                                       vec4(),
                                                       100))); //Inferior
        scene.hit_vector.push_back(new plane(vec4(0, 0, -10, 1), vec4(0, 0, 1, 0),
                                             new phong(vec4(0.7, 0.7, 0.7, 0),
                                                       vec4(),
                                                       0))); //Frontal
        scene.hit_vector.push_back(new plane(vec4(5, 0, 0, 1), vec4(-1, 0, 0, 0),
                                             new phong(vec4(0.1, 0.4, 0.1, 0),
                                                       vec4(),
                                                       10))); //Derecho
        scene.hit_vector.push_back(new plane(vec4(-5, 0, 0, 1), vec4(1, 0, 0, 0),
                                             new phong(vec4(0.4, 0.1, 0.1, 0),
                                                       vec4(0., 0, 0, 0),
                                                       0)));//Izquierdo
        scene.hit_vector.push_back(new plane(vec4(0, 0, 1, 1), vec4(0, 0, -1, 0),
                                             new lambertian(vec4(0.8, 0.4, 0.1, 0)))); //Trasera





        vec4 point;
        std::vector<std::array<float, 3>> pixels(opt.width * opt.height);
        //float color_resolution = 0, max_rgb;
        clock_t start, end;
        start = clock();
        //for (int j = 0; j < cam.height; j++) {
        //    for (int i = 0; i < cam.width; i++) {
        //        rgb = vec4(0, 0, 0, 0);
//
//
        //        rgb = trace(cam, i, j, scene, lights);
        //        if (max(rgb) > color_resolution) {
        //            color_resolution = max_rgb;
        //        }
        //        pixels.push_back({rgb.r(), rgb.g(), rgb.b()});
        //    }
        //    std::cout << "\r" << 100 * j / cam.height << "%";
        //}

        int image_size = opt.width * opt.height;

        volatile std::atomic<std::size_t> count(0), color_resolution(0), max_rgb(0);
        std::vector<std::future<void>> future_vector;
        std::array<float, 3> pixel{};
        vec4 rgb = vec4();
        while (opt.cores--) {
            future_vector.emplace_back(
                    std::async([=, &scene, &count, &color_resolution, &max_rgb, &pixels]() {
                        while (true) {
                            int index = count++;
                            if (index % (image_size / 100) == 0)
                                std::cout << "\r" << 100 * index / image_size << "%";
                            if (index >= image_size) {
                                break;
                            }
                            int x = index % opt.width;
                            int y = index / opt.width;
                            vec4 rgb = trace(cam, x, y, scene, lights);
                            pixels[index] = {rgb.r(), rgb.g(), rgb.b()};;
                            if (max(rgb) > max_rgb) {
                                max_rgb = max(rgb);
                                color_resolution = max_rgb;
                            }
                        }
                    }));
        }
        //std::vector<std::future<void>> future_vector2= future_vector;
        future_vector[0].wait();
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
        image i3 = equalize(gamma(clamp(hdr, 6000), 2), 1023);
        i3.save("image_fix3.ppm");

    } else {
        image i = image("image_hdr.ppm");
        i = gamma(clamp(i, 1000), 2.2);
        i.save("image_fix4.ppm");
    }
    return 0;
}
