//
// Created by muniia on 19/09/19.
//

#include <array>
#include <cmath>

#ifndef TRABAJO1_VEC3_H
#define TRABAJO1_VEC3_H

#endif //TRABAJO1_VEC3_H

/**
 * Clase que define vectores, puntos, etc
 */
class vec4 {
public:
    std::array<float, 4> vector;

    vec4();

    vec4(float x, float y, float z, float w);
    //vec4(const vec4 &vect);

    float x() const;

    float y() const;

    float z() const;

    float w() const;

    float r() const;

    float g() const;

    float b() const;

    //std::array<float, 4> vec() const { return vector; };
};

class mat4 {
    std::array<std::array<float, 4>, 4> matrix;

public:
    mat4(vec4 i, vec4 j, vec4 k, vec4 x) : matrix({i.x(), j.x(), k.x(), x.x(),
                                                   i.y(), j.y(), k.y(), x.y(),
                                                   i.z(), j.z(), k.z(), x.z(),
                                                   i.w(), j.w(), k.w(), x.w()}) {}

    mat4();


    mat4(std::array<std::array<float, 4>, 4> m) : matrix(m) {}

    std::array<std::array<float, 4>, 4> mat() const {
        return matrix;
    }

};

vec4::vec4() : vector({0, 0, 0, 0}) {};

vec4::vec4(float x, float y, float z, float w) : vector({x, y, z, w}) {}

inline float vec4::x() const {
    return vector[0];
}

inline float vec4::y() const {
    return vector[1];
}

inline float vec4::z() const {
    return vector[2];
}

inline float vec4::w() const {
    return vector[3];
}

inline float vec4::r() const {
    return vector[0];
}

inline float vec4::g() const {
    return vector[1];
}

inline float vec4::b() const {
    return vector[2];
}

inline float modulus(vec4 v) {
    return sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}

/**
 * Producto escalar
 * @param d
 * @param e
 * @return
 */
inline float dot(const vec4 d, const vec4 e) {
    return d.x() * e.x() +
           d.y() * e.y() +
           d.z() * e.z();
}

/**
 * Producto vectorial
 * @param d
 * @param e
 * @return
 */
inline vec4 cross(const vec4 d, const vec4 e) {
    return {d.y() * e.z() - d.z() * e.y(),
            d.z() * e.x() - d.x() * e.z(),
            d.x() * e.y() - d.y() * e.x(),
            0};
}


inline vec4 operator+(const vec4 &v1, const vec4 &v2) {
    return {v1.x() + v2.x(),
            v1.y() + v2.y(),
            v1.z() + v2.z(),
            v1.w() + v2.w()};
}

inline vec4 operator-(const vec4 &v1, const vec4 &v2) {
    return {v1.x() - v2.x(),
            v1.y() - v2.y(),
            v1.z() - v2.z(),
            v1.w() - v2.w()};
}

inline vec4 operator*(const vec4 &v1, const vec4 &v2) {
    return {v1.x() * v2.x(),
            v1.y() * v2.y(),
            v1.z() * v2.z(),
            v1.w() * v2.w()};
}

inline vec4 operator*(float s, const vec4 &v) {
    return vec4(s * v.x(),
                s * v.y(),
                s * v.z(),
                v.w());
}

inline vec4 operator*(const vec4 &v, float s) {
    return vec4(s * v.x(),
                s * v.y(),
                s * v.z(),
                v.w());
}

inline vec4 operator/(const vec4 &v, float s) {
    return vec4(v.x() / s,
                v.y() / s,
                v.z() / s,
                v.w());
}

std::ostream &operator<<(std::ostream &os, const vec4 &v) {
    return os << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.w();
}

inline bool operator==(const vec4 &v1, const vec4 &v2) {
    float epsilon = 10e-6;
    return v1.x() - v2.x() < epsilon &&
           v1.y() - v2.y() < epsilon &&
           v1.z() - v2.z() < epsilon &&
           v1.w() - v2.w() < epsilon;
}

inline vec4 operator+=(const vec4 &v1, const vec4 &v2) {
    return v1 + v2;
}

inline vec4 normalize(vec4 v) {
    return v / modulus(v); //vec4(v.x() / modulus, v.y() / modulus, v.z() / modulus);
}

inline vec4 normalize_point(vec4 v) {
    return vec4(v.x() / v.w(),
                v.y() / v.w(),
                v.z() / v.w(),
                1);
}

inline float max(vec4 v) {
    return std::max(v.b(), std::max(v.r(), v.g()));
}

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline mat4 identity() {
    return mat4({1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1});
}

inline mat4 translation(float x, float y, float z) {
    return mat4({1, 0, 0, x,
                 0, 1, 0, y,
                 0, 0, 1, z,
                 0, 0, 0, 1});
}

inline mat4 translation(vec4 vec) {
    return mat4({1, 0, 0, vec.x(),
                 0, 1, 0, vec.y(),
                 0, 0, 1, vec.z(),
                 0, 0, 0, 1});
}

inline mat4 scale(float x, float y, float z) {
    return mat4({x, 0, 0, 0,
                 0, y, 0, 0,
                 0, 0, z, 0,
                 0, 0, 0, 1});
}

inline mat4 rotation_x(float theta) {
    return mat4({1, 0, 0, 0,
                 0, std::cos(theta), -std::sin(theta), 0,
                 0, std::sin(theta), std::cos(theta), 0,
                 0, 0, 0, 1});
}

inline mat4 rotation_y(float theta) {
    return mat4({std::cos(theta), 0, std::sin(theta), 0,
                 0, 1, 0, 0,
                 -std::sin(theta), 0, std::cos(theta), 0,
                 0, 0, 0, 1});
}

inline mat4 rotation_z(float theta) {
    return mat4({std::cos(theta), -std::sin(theta), 0, 0,
                 std::sin(theta), std::cos(theta), 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1});
}

inline mat4 change_base(vec4 u, vec4 v, vec4 w, vec4 o) {
    return mat4({u.x(), v.x(), w.x(), o.x(),
                 u.y(), v.y(), w.y(), o.y(),
                 u.z(), v.z(), w.z(), o.z(),
                 0, 0, 0, 1});
}

/**
 * Multiplicacion matrix por vector
 * @param M
 * @param v
 * @param w
 * @return
 */
inline vec4 matmul(const mat4 &M, const vec4 &v) {
    //std::array<std::array<float, 4>, 4> Mat1 = M.mat();
    //std::array<float, 4> vec = v.vec();
    std::array<float, 4> result;

    for (int i = 0; i < 4; i++) {
        result[i] = 0;
        for (int j = 0; j < 4; j++) {
            result[i] += M.mat()[i][j] * v.vector[j];
        }
    }
    if (v.w() == 1) {
        result[0] /= result[3];
        result[1] /= result[3];
        result[2] /= result[3];
    }

    return vec4(result[0], result[1], result[2], result[3]);
}

/**
 * Multiplicación de matrices
 * @param M1
 * @param M2
 * @return
 */
inline mat4 matmul(const mat4 &M1, const mat4 &M2) {
    //std::array<std::array<float, 4>, 4> mat1 = M1.mat();
    //std::array<std::array<float, 4>, 4> mat2 = M2.mat();
    std::array<std::array<float, 4>, 4> result;

    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            result[i][j] = 0;
            for (int k = 0; k <= 3; k++) {
                result[i][j] = result[i][j] + M1.mat()[i][k] * M2.mat()[k][j];
            }
        }
    }
    return mat4(result);
}

/**
 * Multiplicación de matrices
 * @param M1
 * @param M2
 * @return
 */
inline mat4 operator*(const mat4 &M1, const mat4 &M2) {
    std::array<std::array<float, 4>, 4> Mat1 = M1.mat();
    std::array<std::array<float, 4>, 4> Mat2 = M2.mat();
    std::array<std::array<float, 4>, 4> result;

    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            result[i][j] = 0;
            for (int k = 0; k <= 3; k++) {
                result[i][j] = result[i][j] + Mat1[i][k] * Mat2[k][j];
            }
        }
    }
    return mat4(result);
}

/**
 * Multiplicación matriz vector
 * @param M
 * @param v
 * @return
 */
inline vec4 operator*(const mat4 &M, const vec4 &v) {
    std::array<std::array<float, 4>, 4> Mat1 = M.mat();
    std::array<float, 4> vector = v.vector;
    std::array<float, 4> result;

    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            result[i] = 0;
            for (int k = 0; k <= 3; k++) {
                result[i] += Mat1[i][k] * vector[k];
            }
        }
    }
    return vec4(result[0], result[1], result[2], result[3]);
}

inline std::ostream &operator<<(std::ostream &os, const mat4 &m) {
    std::string mat_string;

    for (int i = 0; i <= 3; i++) {

        for (int j = 0; j <= 3; j++) {
            mat_string += std::to_string(m.mat()[i][j]) + "\t";
        }
        mat_string += "\n";
    }
    return os << mat_string;

}

inline std::array<std::array<float, 3>, 3> inverse3x3(std::array<std::array<float, 3>, 3> M3) {
    std::array<std::array<float, 3>, 3> inv;
    float det = M3[0][0] * (M3[1][1] * M3[2][2] - M3[2][1] * M3[1][2]) -
                M3[0][1] * (M3[1][0] * M3[2][2] - M3[1][2] * M3[2][0]) +
                M3[0][2] * (M3[1][0] * M3[2][1] - M3[1][1] * M3[2][0]);

    float invdet = 1 / det;

    inv[0][0] = (M3[1][1] * M3[2][2] - M3[2][1] * M3[1][2]) * invdet;
    inv[0][1] = (M3[0][2] * M3[2][1] - M3[0][1] * M3[2][2]) * invdet;
    inv[0][2] = (M3[0][1] * M3[1][2] - M3[0][2] * M3[1][1]) * invdet;
    inv[1][0] = (M3[1][2] * M3[2][0] - M3[1][0] * M3[2][2]) * invdet;
    inv[1][1] = (M3[0][0] * M3[2][2] - M3[0][2] * M3[2][0]) * invdet;
    inv[1][2] = (M3[1][0] * M3[0][2] - M3[0][0] * M3[1][2]) * invdet;
    inv[2][0] = (M3[1][0] * M3[2][1] - M3[2][0] * M3[1][1]) * invdet;
    inv[2][1] = (M3[2][0] * M3[0][1] - M3[0][0] * M3[2][1]) * invdet;
    inv[2][2] = (M3[0][0] * M3[1][1] - M3[1][0] * M3[0][1]) * invdet;


    return inv;
}

inline mat4 inverse(mat4 M) {
    //std::array<std::array<float, 4>, 4> inv;
    std::array<std::array<float, 3>, 3> R = {M.mat()[0][0], M.mat()[0][1], M.mat()[0][2],
                                             M.mat()[1][0], M.mat()[1][1], M.mat()[1][2],
                                             M.mat()[2][0], M.mat()[2][1], M.mat()[2][2]};
    std::array<std::array<float, 3>, 3> Rinv = inverse3x3(R);

    return mat4({Rinv[0][0], Rinv[0][1], Rinv[0][2],
                 -(Rinv[0][0] * M.mat()[0][3] + Rinv[0][1] * M.mat()[1][3] + Rinv[0][2] * M.mat()[2][3]),
                 Rinv[1][0], Rinv[1][1], Rinv[1][2],
                 -(Rinv[1][0] * M.mat()[0][3] + Rinv[1][1] * M.mat()[1][3] + Rinv[1][2] * M.mat()[2][3]),
                 Rinv[2][0], Rinv[2][1], Rinv[2][2],
                 -(Rinv[2][0] * M.mat()[0][3] + Rinv[2][1] * M.mat()[1][3] + Rinv[2][2] * M.mat()[2][3]),
                 0, 0, 0, 1});


}