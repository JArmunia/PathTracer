//
// Created by muniia on 21/10/19.
//

#ifndef TRABAJO2_IMAGE_H
#define TRABAJO2_IMAGE_H

#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <array>
#include <iostream>
#include <regex>
#include <tgmath.h>

class image {
public:
    std::string formatIdentification;
    //float MAX;
    int resolution[2];
    float color_resolution;
    std::vector<std::array<float, 3>> pixels;

    image(const std::string &file);

    image(std::string fId, int resolution_x, int resolution_y,
          float color_resolution, std::vector<std::array<float, 3>> pixels);

    void save(const std::string &file);
};

image::image(std::string fId, int resolution_x, int resolution_y, float c_r,
             std::vector<std::array<float, 3>> px) {
    formatIdentification = fId;
    resolution[0] = resolution_x;
    resolution[1] = resolution_y;
    color_resolution = c_r;
    pixels = px;
}

image::image(const std::string &file) {
    std::string line;
    std::ifstream myfile(file);
    std::string value;
    bool finished = false;
    bool hasMAX = false;
    float MAX;

    if (myfile.is_open()) {
        std::regex comment("^#.*");
        std::regex max("^#MAX=.*");
        std::regex number("[0-9]+");
        int i = 0;
        std::smatch match;
        while (!finished) {
            getline(myfile, line);
            if (std::regex_match(line, comment)) {
                if (std::regex_search(line, max)) {
                    std::regex_search(line, match, number);
                    hasMAX = true;
                    MAX = stof(match[0]);
                }
            } else {
                switch (i) {
                    case 0:
                        formatIdentification = line;
                        i++;
                        break;

                    case 1:

                        std::regex_search(line, match, number, std::regex_constants::match_any);
                        resolution[0] = stoi(match[0]);
                        resolution[1] = stoi(match.suffix());
                        i++;
                        break;

                    case 2:
                        std::regex_search(line, match, number);
                        color_resolution = stof(match[0]);
                        i++;
                        finished = true;
                        break;

                    default:

                        break;
                }
            }

        }
        //std::cout << "Format: " << formatIdentification << " MAX: " << MAX << " resolution: " << resolution[0] << ", "
        //          << resolution[1] << std::endl;

        float v;
        for (i = 0; i < resolution[0] * resolution[1]; i++) {
            std::array<float, 3> pixel{};
            for (int c = 0; c < 3; c++) {

                myfile >> v;
                if (hasMAX) pixel[c] = MAX * v / color_resolution;
                else pixel[c] = v;

            }
            pixels.push_back(pixel);
        }
        if (hasMAX)color_resolution = MAX;
    }

    myfile.close();
}

void image::save(const std::string &file) {
    std::ofstream myfile;
    myfile.open(file);
    myfile << formatIdentification << std::endl;
    myfile << resolution[0] << " " << resolution[1] << std::endl;
    myfile << std::to_string(color_resolution) << std::endl;

    for (int j = 0; j < resolution[1]; j++) {
        for (int i = 0; i < resolution[0]; i++) {
            std::array<float, 3> pixel = pixels[j * resolution[0] + i];

            myfile << std::to_string(pixel[0]) << " "
                   << std::to_string(pixel[1]) << " "
                   << std::to_string(pixel[2]) << "     ";
        }
        myfile << std::endl;
    }

    myfile.close();
}


image clamp(image im, int c = 1) {
    for (int j = 0; j < im.resolution[1]; j++) {
        for (int i = 0; i < im.resolution[0]; i++) {
            std::array<float, 3> pixel = im.pixels[j * im.resolution[0] + i];
            for (int color = 0; color < 3; color++) {
                if (pixel[color] >= c) {
                    im.pixels[j * im.resolution[0] + i][color] = c;
                }
            }
        }
    }
    im.color_resolution = c;
    return im;
}

image equalize(image im, int c) {
    for (int j = 0; j < im.resolution[1]; j++) {
        for (int i = 0; i < im.resolution[0]; i++) {
            std::array<float, 3> pixel = im.pixels[j * im.resolution[0] + i];
            for (int color = 0; color < 3; color++) {
                im.pixels[j * im.resolution[0] + i][color] = (pixel[color] / im.color_resolution) * c;
            }
        }
    }

    im.color_resolution = c;
    return im;
}


image gamma(image im, float gamma) {
    for (int j = 0; j < im.resolution[1]; j++) {
        for (int i = 0; i < im.resolution[0]; i++) {
            std::array<float, 3> pixel = im.pixels[j * im.resolution[0] + i];
            for (int color = 0; color < 3; color++) {
                im.pixels[j * im.resolution[0] + i][color] = std::pow(pixel[color], (1 / gamma));
            }
        }
    }
    im.color_resolution = std::pow(im.color_resolution, (1 / gamma));
    return im;
}


#endif //TRABAJO2_IMAGE_H
