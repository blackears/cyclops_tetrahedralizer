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

 
#ifndef CYCLOPS_MATH_H
#define CYCLOPS_MATH_H

#include  <vector>

namespace CyclopsTetra3D {

template<typename T = float>
struct Vector3 {
    T x;
    T y;
    T z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    T magnitude_squared() const { return x * x + y * y + z * z; }
    T magnitude() const { return sqrt(x * x + y * y + z * z); }
    Vector3<T> cross(const Vector3<T> rhs) const { return Vector3<T>(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x); }
    T dot(const Vector3<T> rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }

    Vector3<T> min(const Vector3<T> rhs) const { return Vector3<T>(min(x, rhs.x), min(y, rhs.y), min(z, rhs.z)); }
    Vector3<T> max(const Vector3<T> rhs) const { return Vector3<T>(max(x, rhs.x), max(y, rhs.y), max(z, rhs.z)); }

    Vector3<T>& operator+=(const Vector3<T>& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
    }

    friend Vector3<T> operator+(Vector3<T> lhs, const Vector3<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    Vector3<T>& operator-=(const Vector3<T>& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        return *this;
    }

    friend Vector3<T> operator-(Vector3<T> lhs, const Vector3<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Vector3<T>& operator*=(T rhs) {
        this->x *= rhs;
        this->y *= rhs;
        this->z *= rhs;
        return *this;
    }

    friend Vector3<T> operator*(Vector3<T> lhs, T rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3<T>& operator/=(T rhs) {
        this->x /= rhs;
        this->y /= rhs;
        this->z /= rhs;
        return *this;
    }

    friend Vector3<T> operator/(Vector3<T> lhs, T rhs) {
        lhs /= rhs;
        return lhs;
    }


    Vector3<T>& operator*=(Vector3<T> rhs) {
        this->x *= rhs.x;
        this->y *= rhs.y;
        this->z *= rhs.z;
        return *this;
    }

    friend Vector3<T> operator*(Vector3<T> lhs, Vector3<T> rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3<T>& operator/=(Vector3<T> rhs) {
        this->x /= rhs.x;
        this->y /= rhs.y;
        this->z /= rhs.z;
        return *this;
    }

    friend Vector3<T> operator/(Vector3<T> lhs, Vector3<T> rhs) {
        lhs /= rhs;
        return lhs;
    }

};


template<typename T>
struct BoundingBox {
    Vector3<T> bb_min;
    Vector3<T> bb_max;

    BoundingBox() : bb_min(Vector3<T>()), bb_max(Vector3<T>()) {}
    BoundingBox(Vector3<T> bb_min, Vector3<T> bb_max) : bb_min(bb_min), bb_max(bb_max) {}

    BoundingBox<T> merge(const BoundingBox<T>& other) const {
        Vector3<T> new_bb_min = Vector3<T>(std::min(bb_min.x, other.bb_min.x),
            std::min(bb_min.y, other.bb_min.y),
            std::min(bb_min.z, other.bb_min.z));
        Vector3<T> new_bb_max = Vector3<T>(std::max(bb_max.x, other.bb_max.x),
            std::max(bb_max.y, other.bb_max.y),
            std::max(bb_max.z, other.bb_max.z));
        return BoundingBox<T>(new_bb_min, new_bb_max);
    }

    Vector3<T> center() const {
        return (bb_min + bb_max) / 2.0;
    }

    Vector3<T> size() const {
        return bb_max - bb_min;
    }

    bool intersects_ray(Vector3<T> ray_origin, Vector3<T> ray_direction) {
        //Slab method
        //https://en.wikipedia.org/wiki/Slab_method
        Vector3<T> t_low = (bb_min - ray_origin) / ray_direction;
        Vector3<T> t_high = (bb_max - ray_origin) / ray_direction;

        Vector3<T> t_close = t_low.min(t_high);
        Vector3<T> t_far = t_low.max(t_high);

        T t_close_max = std::max(std::max(t_close.x, t_close.y), t_close.z);
        T t_far_min = std::min(std::min(t_far.x, t_far.y), t_far.z);
        return t_close_max <= t_far_min;
    }
};


}

#endif //CYCLOPS_MATH_H