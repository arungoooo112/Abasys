#ifndef TINYIGA_COORD_H
#define TINYIGA_COORD_H

#include "glm/glm.hpp"

template<typename T>
using vec2 = glm::vec<2, T>;

template<typename T>
using vec3 = glm::vec<3, T>;

template<typename T>
using hvec3 = glm::vec<4, T>;

using vec3d = vec3<double>;
using hvec3d = hvec3<double>;

using idx2 = vec2<int>;

template<typename T>
vec3<T> rotateZ(const vec3<T>& a, T alpha) {
    T tx = a.x;
    T ty = a.y;
    T x = tx * cos(alpha) - ty * sin(alpha);
    T y = tx * sin(alpha) + ty * cos(alpha);
    return vec3<T>(x, y, a.z);
}

#endif //TINYIGA_COORD_H