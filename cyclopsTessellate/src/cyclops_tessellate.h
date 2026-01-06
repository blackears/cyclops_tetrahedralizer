/*
 * Copyright (c) 2026 Mark McKay
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef CYCLOPS_TESSELLATE_H
#define CYCLOPS_TESSELLATE_H

namespace CyclopsTessellate3D {

//typedef float real_t;

template<typename T = float>
struct Vector3 {
    T x;
    T y;
    T z;

    Vector3<T> cross(const Vector3<T> a) const { return Vector3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
    T dot(const Vector3<T> a) const { return x * a.x + y * a.y + z * a.z; }

    Vector3<T>& operator+=(const Vector3<T>& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }

    friend Vector3<T> operator+(Vector3<T> lhs, const Vector3<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    Vector3<T>& operator-=(const Vector3<T>& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
    }

    friend Vector3<T> operator-(Vector3<T> lhs, const Vector3<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }
};

struct Triangle {
    int v0_idx;
    int v1_idx;
    int v2_idx;
};

template<typename T = float>
struct Tetrahedron {
    int v0_idx;
    int v1_idx;
    int v2_idx;
    int v3_idx;

    Vector3<T> calc_circum_center(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2, Vector3<T> v3) const;
};

// class Point3D {
//     real_t x;
//     real_t y;
//     real_t z;
// };



template<typename T = float>
class CyclopsTess3D {
    Vector3<T>* point_list;

public:

    void tessellate_tetrahedra(float* mesh_points);

};

}

#endif //CYCLOPS_TESSELLATE_H