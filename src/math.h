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

typedef float real;

struct Vector3 {
    real x;
    real y;
    real z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(real x, real y, real z) : x(x), y(y), z(z) {}

    real magnitude_squared() const { return x * x + y * y + z * z; }
    real magnitude() const { return sqrt(x * x + y * y + z * z); }
    Vector3 cross(const Vector3 rhs) const { return Vector3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x); }
    real dot(const Vector3 rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
    Vector3 normalized() const {
        real mag = magnitude();
        if (mag == 0) return Vector3(0, 0, 0);
        return Vector3(x / mag, y / mag, z / mag);
    }

    Vector3 min(const Vector3 rhs) const { return Vector3(std::min(x, rhs.x), std::min(y, rhs.y), std::min(z, rhs.z)); }
    Vector3 max(const Vector3 rhs) const { return Vector3(std::max(x, rhs.x), std::max(y, rhs.y), std::max(z, rhs.z)); }

    real& operator[](int index) {
        if (index == 0) return x;
        else if (index == 1) return y;
        else return z;
    }

    Vector3& operator+=(const Vector3& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
    }

    friend Vector3 operator+(Vector3 lhs, const Vector3& rhs) {
        lhs += rhs;
        return lhs;
    }

    Vector3& operator-=(const Vector3& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        return *this;
    }

    friend Vector3 operator-(Vector3 lhs, const Vector3& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Vector3& operator*=(real rhs) {
        this->x *= rhs;
        this->y *= rhs;
        this->z *= rhs;
        return *this;
    }

    friend Vector3 operator*(Vector3 lhs, real rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3& operator/=(real rhs) {
        this->x /= rhs;
        this->y /= rhs;
        this->z /= rhs;
        return *this;
    }

    friend Vector3 operator/(Vector3 lhs, real rhs) {
        lhs /= rhs;
        return lhs;
    }


    Vector3& operator*=(Vector3 rhs) {
        this->x *= rhs.x;
        this->y *= rhs.y;
        this->z *= rhs.z;
        return *this;
    }

    friend Vector3 operator*(Vector3 lhs, Vector3 rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3& operator/=(Vector3 rhs) {
        this->x /= rhs.x;
        this->y /= rhs.y;
        this->z /= rhs.z;
        return *this;
    }

    friend Vector3 operator/(Vector3 lhs, Vector3 rhs) {
        lhs /= rhs;
        return lhs;
    }

};


struct BoundingBox {
    Vector3 bb_min;
    Vector3 bb_max;

    BoundingBox() : bb_min(Vector3()), bb_max(Vector3()) {}
    BoundingBox(Vector3 bb_min, Vector3 bb_max) : bb_min(bb_min), bb_max(bb_max) {}
    
    BoundingBox merge(const BoundingBox& other) const {
        Vector3 new_bb_min = Vector3(std::min(bb_min.x, other.bb_min.x),
            std::min(bb_min.y, other.bb_min.y),
            std::min(bb_min.z, other.bb_min.z));
        Vector3 new_bb_max = Vector3(std::max(bb_max.x, other.bb_max.x),
            std::max(bb_max.y, other.bb_max.y),
            std::max(bb_max.z, other.bb_max.z));
        return BoundingBox(new_bb_min, new_bb_max);
    }

    Vector3 center() const {
        return (bb_min + bb_max) / 2.0;
    }

    Vector3 size() const {
        return bb_max - bb_min;
    }

    bool intersects_ray(const Vector3& ray_origin, const Vector3& ray_direction) const {
        //Slab method
        //https://en.wikipedia.org/wiki/Slab_method
        Vector3 t_low = (bb_min - ray_origin) / ray_direction;
        Vector3 t_high = (bb_max - ray_origin) / ray_direction;
        Vector3 t_close = t_low.min(t_high);
        Vector3 t_far = t_low.max(t_high);

        real t_close_max = std::max(std::max(t_close.x, t_close.y), t_close.z);
        real t_far_min = std::min(std::min(t_far.x, t_far.y), t_far.z);
        return t_close_max <= t_far_min;
    }
};


struct Plane {
    Vector3 normal;
    //distance along normal from origin to plane
    real dot_origin;

    Plane() : normal(Vector3()), dot_origin(0) {}
    Plane(const Vector3& normal, const Vector3& origin) : normal(normal), dot_origin(origin.dot(normal)) {}
    Plane(const Vector3& p0, const Vector3& p1, const Vector3& p2) : normal((p1 - p0).cross(p2 - p0)), dot_origin(p0.dot(normal)) {}
    
    real distance_to_plane(const Vector3& p) const {
        return normal.dot(p) - dot_origin;
    }

    bool intersect_ray(const Vector3& ray_origin, const Vector3& ray_direction, Vector3& out_intersection) const {
        real denom = normal.dot(ray_direction);
        if (denom == 0.0) {
            return false;
        }
        real numer = dot_origin - normal.dot(ray_origin);

        real s = numer / denom;
        out_intersection = ray_origin + ray_direction * s;
        return true;
    }
};

}

#endif //CYCLOPS_MATH_H