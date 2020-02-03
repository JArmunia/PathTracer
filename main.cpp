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
#include "rectangle.h"


struct options {
    int width = 450;
    int height = 500;
    float h_fov = 120;

    int rays_per_pixel = 20;
    float shadow_bias = 10e-4;
    float intersection_bias = 10e-4;
    int cores = std::thread::hardware_concurrency();
    int scene = 6;


} opt;

vec4 color(ray r, hitable_list world, const std::vector<point_light *> &lights, int n) {
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
                vec4 light_direction = normalize(light->position - rec.p);
                if (event == SPECULAR) {
                    vec4 a();
                }
                rgb = rgb + BRDF(event, rec, r.direction, light_direction,
                                 luminance(*light, rec.p)) *
                            std::abs(dot(normalize(rec.normal), normalize(light_direction)));
            } else {
                // rgb += 0
            }
        }
        if (rec.mat->isLight)
            return rec.mat->getLight(rec);
        vec4 next_direction;
        vec4 point;
        if (event != ABSORTION) {
            switch (event) {
                case DIFFUSE: {
                    next_direction = cosine_sampling_random_direction(rec);
                    //rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, next_direction);
                    break;
                }

                case PHONG_SPECULAR: {
                    //scattered = ray(rec.p, cosine_sampling_random_direction(rec));
                    next_direction = cosine_sampling_random_direction(rec);
                    //next_direction = lobe_sampling_random_direction(rec, reflected(r.direction, rec.normal));
                    //rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, lobe_sampling_random_direction(rec));
                    break;
                }

                case SPECULAR: {
                    //vec4 next_direction = r.next_direction - 2 * rec.normal * dot(r.next_direction, rec.normal);
                    //scattered = ray(rec.p, reflected(r.next_direction, rec.normal));
                    next_direction = reflected(r.direction, rec.normal);
                    //rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    break;
                }

                case REFRACTION: {
                    next_direction = refract(r.direction, rec.normal, rec.mat->refraction_index);
                    //rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
                    //scattered = ray(rec.p, next_direction);
                    break;
                }
            }
            rec.p = rec.p + normalize(next_direction) * opt.intersection_bias;
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


        hitable_list scene;

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// SPHERES ///////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        vec4 white = vec4(.85, .85, .85, 0);
        vec4 red = vec4(.85, .085, .085, 0);
        vec4 green = vec4(.085, .85, .085, 0);
        vec4 orange = vec4(.85, .6, .02, 0);
        vec4 blue = vec4(.085, .085, .85, 0);

        sphere pixar = sphere(vec4(0, 0, 0, 1), 1.5,
                              new texture(image("pixar_ball.ppm"), 2.f, 1));
        pixar.basis = pixar.basis * translation(vec4(3.5, -3.5, -8.5, 0)) * rotation_x(-M_PI / 5) *
                      rotation_y(M_PI / 3);

        switch (opt.scene) {
            case 1:
                //lights.push_back(new point_light(vec4(0, 2, -4, 1),
                //                                 vec4(1, 1, 1, 0), 60000));

                scene.hit_vector.push_back(new sphere(vec4(-3, -1, -6, 1), 1.,
                                                      new glass(vec4(0.9, 0.9, 0.9, 0), 1.8)));
//
                scene.hit_vector.push_back(new sphere(vec4(0, -2.5, -6, 1), 1.5,
                        //new specular(vec4(0.9, 0.9, 0.9, 0))));
                                                      new phong(vec4(0.5, 0.5, 0.5, 0), vec4(), 3)));

                scene.hit_vector.push_back(new sphere(vec4(2, 2.5, -6, 1), 1.,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));
                scene.hit_vector.push_back(new sphere(vec4(-2, 2.5, -6, 1), 1.,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));

                scene.hit_vector.push_back(new sphere(vec4(3, -2.5, -8, 1), 1.5,
                                                      new phong(vec4(0.1, 0.6, 0.6, 0),
                                                                vec4(0., 0, 0., 0),
                                                                15)));


                scene.hit_vector.push_back(new plane(vec4(0, 4, 0, 1), vec4(0, -1, 0, 0),
                        //new light(vec4(1, 1, 1, 0), 10000))); //Superior
                                                     new phong(vec4(0.4, 0.4, 0.4, 0), vec4(), 0)));

                scene.hit_vector.push_back(new plane(vec4(0, -4, 0, 1), vec4(0, 1, 0, 0),
                                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                                                               vec4(),
                                                               100))); //Inferior
                //scene.hit_vector.push_back(new plane(vec4(0, 0, -10, 1), vec4(0, 0, 1, 0),
                //                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                //                                               vec4(),
                //                                               0))); //Frontal

                scene.hit_vector.push_back(new finite_plane(vec4(0, 0, -10, 1), 10, 8,
                                                            new texture(image("../Textures/moon.ppm"), 6000)));
                //new specular(vec4(.9, .9, .9, 0))));
                scene.hit_vector.push_back(new plane(vec4(5, 0, 0, 1), vec4(-1, 0, 0, 0),
                        //new phong(vec4(0.1, 0.4, 0.1, 0),
                        //          vec4(),
                        //          10))); //Derecho
                                                     new specular(vec4(.9, .9, .9, 0))));
                scene.hit_vector.push_back(new plane(vec4(-5, 0, 0, 1), vec4(1, 0, 0, 0),
                                                     new phong(vec4(0.4, 0.1, 0.1, 0),
                                                               vec4(0., 0, 0, 0),
                                                               0)));//Izquierdo
                break;
            case 2:
                // Simple cornell box
                lights.push_back(new point_light(vec4(0, 4.5, -7.5, 1),
                                                 vec4(1, 1, 1, 0), 60000));

                scene.hit_vector.push_back(new plane(vec4(0, 6, 0, 1), vec4(0, -1, 0, 0),
                        //new light(vec4(1, 1, 1, 0), 10000))); //Superior
                                                     new phong(vec4(0.4, 0.4, 0.4, 0), vec4(), 0)));

                scene.hit_vector.push_back(new plane(vec4(0, -6, 0, 1), vec4(0, 1, 0, 0),
                                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                                                               vec4(),
                                                               100))); //Inferior

                scene.hit_vector.push_back(new plane(vec4(0, 0, -12, 1), vec4(0, 0, 1, 0),
                                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                                                               vec4(),
                                                               0))); //Frontal

                scene.hit_vector.push_back(new plane(vec4(6, 0, 0, 1), vec4(-1, 0, 0, 0),
                                                     new phong(vec4(0.1, 0.7, 0.1, 0),
                                                               vec4(),
                                                               10))); //Derecho

                scene.hit_vector.push_back(new plane(vec4(-6, 0, 0, 1), vec4(1, 0, 0, 0),
                                                     new phong(vec4(0.7, 0.1, 0.1, 0),
                                                               vec4(0., 0, 0, 0),
                                                               0)));//Izquierdo

                scene.hit_vector.push_back(new plane(vec4(0, 0, 1, 1), vec4(0, 0, -1, 0),
                                                     new phong(vec4(0.4, 0.4, 0.4, 0),
                                                               vec4(0., 0, 0, 0),
                                                               0))); //Trasera

                //scene.hit_vector.push_back(new sphere(vec4(3.5, -4, -10, 1), 2,
                //                                      new phong(vec4(0.5, 0.2, 0.2, 0),
                //                                                vec4(0.2, 0.2, 0.2, 0), 10)));

                scene.hit_vector.push_back(new sphere(vec4(3.5, -4, -10, 1), 2,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));

                scene.hit_vector.push_back(new sphere(vec4(-1, -4, -7, 1), 1.5,
                                                      new glass(vec4(0.9, 0.9, 0.9, 0),
                                                                1.8)));
                //pixar.mat->Ks = vec4(0.2,0.2,0.2,0);
                //pixar.mat->alpha = 10;

                scene.hit_vector.push_back(&pixar);

                break;
            case 3 :


                scene.hit_vector.push_back(new plane(vec4(0, 6, 0, 1), vec4(0, -1, 0, 0),
                                                     new light(vec4(1, 1, 1, 0), 10000))); //Superior
                //new phong(vec4(0.4, 0.4, 0.4, 0), vec4(), 0)));

                scene.hit_vector.push_back(new plane(vec4(0, -6, 0, 1), vec4(0, 1, 0, 0),
                                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                                                               vec4(),
                                                               100))); //Inferior

                scene.hit_vector.push_back(new plane(vec4(0, 0, -12, 1), vec4(0, 0, 1, 0),
                                                     new phong(vec4(0.7, 0.7, 0.7, 0),
                                                               vec4(),
                                                               0))); //Frontal

                scene.hit_vector.push_back(new plane(vec4(6, 0, 0, 1), vec4(-1, 0, 0, 0),
                                                     new phong(vec4(0.1, 0.7, 0.1, 0),
                                                               vec4(),
                                                               10))); //Derecho

                scene.hit_vector.push_back(new plane(vec4(-6, 0, 0, 1), vec4(1, 0, 0, 0),
                                                     new phong(vec4(0.7, 0.1, 0.1, 0),
                                                               vec4(0., 0, 0, 0),
                                                               0)));//Izquierdo

                scene.hit_vector.push_back(new plane(vec4(0, 0, 1, 1), vec4(0, 0, -1, 0),
                                                     new phong(vec4(0.4, 0.4, 0.4, 0),
                                                               vec4(0., 0, 0, 0),
                                                               0))); //Trasera

                scene.hit_vector.push_back(new sphere(vec4(3.5, -4, -10, 1), 2,
                                                      new phong(vec4(0.5, 0.2, 0.2, 0),
                                                                vec4(0.2, 0.2, 0.2, 0), 10)));

                scene.hit_vector.push_back(new sphere(vec4(-3.5, -4, -10, 1), 2,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));

                scene.hit_vector.push_back(new sphere(vec4(-1, -4, -7, 1), 1.5,
                                                      new glass(vec4(0.9, 0.9, 0.9, 0),
                                                                1.8)));


                break;

            case 4 :
                //BSDF* white = new Lambertian(w, Vector3(.85, .85, .85));
                //BSDF* red = new Lambertian(w, Vector3(.85, .085, .085));
                //BSDF* green = new Lambertian(w, Vector3(.085, .85, .085));
                //BSDF* orange = new Lambertian(w, Vector3(.85, .6, .02));
                //BSDF* blue = new Lambertian(w, Vector3(.085, .085, .85));

                opt.h_fov = 90;

                lights.push_back(new point_light(vec4(0, 0.9, -3, 1),
                                                 vec4(1, 1, 1, 0), 500000));

                scene.hit_vector.push_back(new plane(vec4(0, 1, 0, 1), vec4(0, -1, 0, 0),
                        //new light(white, 10000))); //Superior
                                                     new phong(white, vec4(), 0)));
                scene.hit_vector.push_back(new plane(vec4(0, -1, 0, 1), vec4(0, 1, 0, 0),
                                                     new phong(white,
                                                               vec4(),
                                                               100))); //Inferior

                scene.hit_vector.push_back(new plane(vec4(0, 0, -4, 1), vec4(0, 0, 1, 0),
                                                     new phong(white,
                                                               vec4(),
                                                               0))); //Frontal

                scene.hit_vector.push_back(new plane(vec4(1, 0, 0, 1), vec4(-1, 0, 0, 0),
                                                     new phong(green,
                                                               vec4(),
                                                               10))); //Derecho

                scene.hit_vector.push_back(new plane(vec4(-1, 0, 0, 1), vec4(1, 0, 0, 0),
                                                     new phong(red,
                                                               vec4(0., 0, 0, 0),
                                                               0)));//Izquierdo


                scene.hit_vector.push_back(new sphere(vec4(0.5, -0.7, -2.5, 1), 0.3,
                                                      new phong(white,
                                                                vec4(0., 0, 0, 0),
                                                                0)));
                scene.hit_vector.push_back(new sphere(vec4(-0.5, -0.5, -1.5, 1), 0.3,
                                                      new phong(red,
                                                                vec4(0., 0, 0, 0),
                                                                0)));
                scene.hit_vector.push_back(new sphere(vec4(0, -0.7, -3, 1), 0.3,
                                                      new phong(white,
                                                                vec4(0., 0, 0, 0),
                                                                0)));

                break;

            case 5 :
                //BSDF* white = new Lambertian(w, Vector3(.85, .85, .85));
                //BSDF* red = new Lambertian(w, Vector3(.85, .085, .085));
                //BSDF* green = new Lambertian(w, Vector3(.085, .85, .085));
                //BSDF* orange = new Lambertian(w, Vector3(.85, .6, .02));
                //BSDF* blue = new Lambertian(w, Vector3(.085, .085, .85));


                lights.push_back(new point_light(vec4(0, 0.9, -3, 1),
                                                 vec4(1, 1, 1, 0), 3000));

                scene.hit_vector.push_back(new plane(vec4(0, 1, 0, 1), vec4(0, -1, 0, 0),
                        //new light(white, 10000))); //Superior
                                                     new phong(white, vec4(), 0)));
                scene.hit_vector.push_back(new plane(vec4(0, -1, 0, 1), vec4(0, 1, 0, 0),
                                                     new phong(white,
                                                               vec4(),
                                                               100))); //Inferior

                scene.hit_vector.push_back(new plane(vec4(0, 0, -4, 1), vec4(0, 0, 1, 0),
                                                     new phong(white,
                                                               vec4(),
                                                               0))); //Frontal

                scene.hit_vector.push_back(new plane(vec4(1, 0, 0, 1), vec4(-1, 0, 0, 0),
                                                     new phong(green,
                                                               vec4(),
                                                               10))); //Derecho

                scene.hit_vector.push_back(new plane(vec4(-1, 0, 0, 1), vec4(1, 0, 0, 0),
                                                     new phong(red,
                                                               vec4(0., 0, 0, 0),
                                                               0)));//Izquierdo


                scene.hit_vector.push_back(new sphere(vec4(0.5, -0.7, -2.5, 1), 0.3,
                                                      new glass(vec4(0.9, 0.9, 0.9, 0), 1.8)));

                scene.hit_vector.push_back(new sphere(vec4(-0.5, -0.5, -1.5, 1), 0.3,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));

                scene.hit_vector.push_back(new sphere(vec4(0, -0.7, -3, 1), 0.3,
                                                      new phong(vec4(0.6, 0.2, 0.2, 0),
                                                                vec4(0.2, 0.2, 0.2, 0),
                                                                10)));

                break;

            case 6:

                lights.push_back(new point_light(vec4(0,0,0,1), vec4(1,1,1,1), 1000));
                finite_plane front(vec4(0, 0, -12, 1), 12, 12,
                        new texture(image("../Textures/Ca_ship_front.ppm")));
                scene.hit_vector.push_back(&front);
                finite_plane top(vec4(0, 6, -6, 1), 12, 12,
                        //new texture(image("../Textures/Ca_ship_top.ppm")));
                        new light(vec4(1,1,1,0), 100));
                top.basis = top.basis * rotation_x(M_PI/2);
                scene.hit_vector.push_back(&top);

                finite_plane left(vec4(-6, 0, -12, 1), 12, 12,
                        new texture(image("../Textures/Ca_ship_front2.ppm")));
                left.basis = left.basis * rotation_y(M_PI / 2);
                scene.hit_vector.push_back(&left);

                finite_plane right(vec4(6, 0, -12, 1), 12, 12,
                        new texture(image("../Textures/Ca_ship_front2.ppm")));
                right.basis = right.basis * rotation_y(-M_PI / 2);
                scene.hit_vector.push_back(&right);

                finite_plane floor(vec4(0, -6, -12, 1), 12, 12,
                        new texture(image("../Textures/lava.ppm"), 5, 5, 50));
                floor.basis = floor.basis * rotation_x(-M_PI/2);
                scene.hit_vector.push_back(&floor);

                //scene.hit_vector.push_back();



                //scene.hit_vector.push_back(new sphere(vec4(3.5, -4, -10, 1), 2,
                //                                      new phong(vec4(0.5, 0.2, 0.2, 0),
                //                                                vec4(0.2, 0.2, 0.2, 0), 10)));

                scene.hit_vector.push_back(new sphere(vec4(-3.5, -4, -10, 1), 2,
                                                      new specular(vec4(0.9, 0.9, 0.9, 0))));

                scene.hit_vector.push_back(new sphere(vec4(-1, -4, -7, 1), 1.5,
                                                      new glass(vec4(0.9, 0.9, 0.9, 0),
                                                                1.8)));
                //pixar.mat->Ks = vec4(0.2,0.2,0.2,0);
                //pixar.mat->alpha = 10;

                scene.hit_vector.push_back(&pixar);

                break;

                /*case 7:
                    //sphere pixar = sphere(vec4(0, 0, 0, 1), 1.5,
                    //                      new texture(vec4(0.9, 0, 0, 0), image("pixar_ball.ppm"), 2.f, 1, 6000));
                    //pixar.basis = pixar.basis * translation(vec4(-3.5, -2.5, -8.5, 0)) * rotation_x(-M_PI / 5) *
                    //              rotation_y(-M_PI / 3);//* rotation_x(-M_PI/2)  ;
                    //scene.hit_vector.push_back(&pixar);


                    scene.hit_vector.push_back(new xy_rect(30, 60, -30, 90, -120,
                                                           new texture(image("building.ppm"), 1.f, 2)));
                    //scene.hit_vector.push_back(new yz_rect(-30, 90, -150, -120, 30,
                    //                                       new texture(image("building.ppm"), 1.f, 2)));

                    finite_plane right_wall(vec4(30, 30, -135, 1), 30, 120,
                                            new texture(image("building.ppm"), 1, 1));
                    right_wall.basis = right_wall.basis * rotation_y(-M_PI / 2);
                    scene.hit_vector.push_back(&right_wall);
                    scene.hit_vector.push_back(new yz_rect(-30, 90, -150, -120, 60,
                                                           new texture(image("building.ppm"), 1.f, 2)));
                    scene.hit_vector.push_back(new xy_rect(30, 60, -30, 90, -150,
                                                           new texture(image("building.ppm"), 1.f, 2)));
    ////
                    scene.hit_vector.push_back(new xz_rect(-300, 300, -300, 300, -30,
                                                           new texture(image("grass.ppm"), 10.f, 20)));
    //
                    scene.hit_vector.push_back(
                            new xz_rect(30, 60, -150, -120, 90, new phong(vec4(0.2, 0.2, 0.2, 0), vec4(), 1)));
                    sphere sky(vec4(0, -0, 0, 1), 5000,
                               new texture(image("../Textures/sunset24.ppm"), 1.f, 2,
                                           90));
                    //new light(vec4(1,1,1,1), 50));
                    sky.basis = sky.basis * rotation_x(-M_PI / 35) * rotation_y(M_PI / 2);
                    scene.hit_vector.push_back(&sky);
                    lights.push_back(new point_light(vec4(400, 100, -450, 1), vec4(1, 1, 1, 1), 20000000));

                    //new phong(vec4(.2,0.3,.9,0),
                    //        vec4(),0)));

                    finite_plane front(vec4(-25, 0, -40, 1), 20, 60,
                                       new texture(image("../Textures/building_dks-1.ppm")));
                    scene.hit_vector.push_back(&front);

                    finite_plane road(vec4(5, -29.9, -250, 1), 30, 500,
                                      new texture(image("../Textures/road2.ppm"), 1, 17, 20));
                    road.basis = road.basis * rotation_x(-M_PI / 2);
                    scene.hit_vector.push_back(&road);
                    break;*/
        }


        //finite_plane fp = finite_plane(vec4(0, 0, -20, 1), 20, 40,
        //                               new texture(vec4(0.9, 0.2, 0.3, 0), image("earthmap.ppm")));
        //fp.basis = fp.basis * rotation_y(-M_PI / 4);
        //scene.hit_vector.push_back(&fp);
        //new texture(vec4(0.9,0.9,0.9,0), image("sky.ppm"))));

        //scene.hit_vector.push_back(new sphere(vec4(0, 5, -6, 1), 1.5,
        //                                      new light(vec4(1, 1, 1, 0), 10000)));
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// PLANES ////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////


        //new specular(vec4(.9, .9, .9, 0))));
        //scene.hit_vector.push_back(new plane(vec4(0, 0, 1, 1), vec4(0, 0, -1, 0),
        //                                     new phong(vec4(0.8, 0.4, 0.1, 0),
        //                                               vec4(0., 0, 0, 0),
        //                                               0))); //Trasera
        //  new light(vec4(1, 1, 1, 0), 10000)));

        //scene.hit_vector.push_back(new box(vec4(-2, -5, -5, 1), vec4(-1, -1, -3, 1),
        //                                   new glass(
        //                                           vec4(0.9, .9, .9, 0), 1.8
        //                                   )));
        //new phong(vec4(0.4, 0.1, 0.1, 0),
        //          vec4(0., 0, 0, 0),
        //          0)));

        //scene.hit_vector.push_back(new box(vec4(-1.9,-4,-4,1), vec4(-1.1,-2,-3.5,1),
        //        //new glass(
        //        //vec4(0.9,.9,.9, 0), 1.8
        //        //)));
        //                                   new phong(vec4(0.4, 0.1, 0.1, 0),
        //                                             vec4(0., 0, 0, 0),
        //                                             0)));






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
