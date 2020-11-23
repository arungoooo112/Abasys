
// 坐标体系主要依赖于glm库，并对其进行了进一步的封装，提供公共的接口。
// 因此你也可以简单的自定义自己的坐标类，
// 需要特别注意的是，自定义的Coord 类，应该能够进行简单的算术运算，如加减，
// 还应该提供一个绕z轴旋转的接口函数
// 坐标值的访问可以满足如下操作：
// vec<double> v(1,2,3);
// cout << v.x << v[0]; //输出坐标x
// cout << v.y << v[1]; //输出坐标x
// cout << v.z << v[2]; //输出坐标x


#ifndef ABASYS_COORD_H
#define ABASYS_COORD_H

#include "glm/glm.hpp"

template<typename T>
using vec3 = glm::vec<3, T>;

template<typename T>
using hvec3 = glm::vec<4, T>;

using vec3f = vec3<float>;
using vec3d = vec3<double>;
using hvec3f = hvec3<float>;
using hvec3d = hvec3<double>;

template<typename T>
vec3<T> rotateZ(const vec3<T>& a, T alpha) {
    T tx = a.x;
    T ty = a.y;
    T x = tx * cos(alpha) - ty * sin(alpha);
    T y = tx * sin(alpha) + ty * cos(alpha);
    return vec3<T>(x, y, a.z);
}

#endif //ABASYS_COORD_H