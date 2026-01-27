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


struct Vector2 {
    real x;
    real y;

    Vector2() : x(0), y(0) {}
    Vector2(real x, real y) : x(x), y(y) {}

    real magnitude_squared() const { return x * x + y * y; }
    real magnitude() const { return sqrt(x * x + y * y); }
    real dot(const Vector2 rhs) const { return x * rhs.x + y * rhs.y; }
    Vector2 normalized() const {
        real mag = magnitude();
        if (mag == 0) return Vector2(0, 0);
        return Vector2(x / mag, y / mag);
    }

    Vector2 min(const Vector2 rhs) const { return Vector2(std::min(x, rhs.x), std::min(y, rhs.y)); }
    Vector2 max(const Vector2 rhs) const { return Vector2(std::max(x, rhs.x), std::max(y, rhs.y)); }

    Vector2 rot_CCW_90() { return Vector2(-y, x); }
    Vector2 rot_CW_90() { return Vector2(y, -x); }

    real operator[](int index) const {
        if (index == 0) return x;
        return y;
    }

    real& operator[](int index) {
        if (index == 0) return x;
        return y;
    }

    Vector2& operator+=(const Vector2& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }

    friend Vector2 operator+(Vector2 lhs, const Vector2& rhs) {
        lhs += rhs;
        return lhs;
    }

    Vector2& operator-=(const Vector2& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
    }

    friend Vector2 operator-(Vector2 lhs, const Vector2& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Vector2& operator*=(real rhs) {
        this->x *= rhs;
        this->y *= rhs;
        return *this;
    }

    friend Vector2 operator*(Vector2 lhs, real rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector2& operator/=(real rhs) {
        this->x /= rhs;
        this->y /= rhs;
        return *this;
    }

    friend Vector2 operator/(Vector2 lhs, real rhs) {
        lhs /= rhs;
        return lhs;
    }


    Vector2& operator*=(Vector2 rhs) {
        this->x *= rhs.x;
        this->y *= rhs.y;
        return *this;
    }

    friend Vector2 operator*(Vector2 lhs, Vector2 rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector2& operator/=(Vector2 rhs) {
        this->x /= rhs.x;
        this->y /= rhs.y;
        return *this;
    }

    friend Vector2 operator/(Vector2 lhs, Vector2 rhs) {
        lhs /= rhs;
        return lhs;
    }

};

struct Vector3 {
    real x;
    real y;
    real z;

    static const Vector3 X_POS;
    static const Vector3 X_NEG;
    static const Vector3 Y_POS;
    static const Vector3 Y_NEG;
    static const Vector3 Z_POS;
    static const Vector3 Z_NEG;

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

    int max_axis() const {
        if (abs(x) > abs(y) && abs(x) > abs(z))
            return 0;
        if (abs(y) > abs(z))
            return 1;
        return 2;
    }

    int min_axis() const {
        if (abs(x) < abs(y) && abs(x) < abs(z))
            return 0;
        if (abs(y) < abs(z))
            return 1;
        return 2;
    }

    real operator[](int index) const {
        if (index == 0) return x;
        else if (index == 1) return y;
        else return z;
    }

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


struct Rectangle {
    Vector2 bb_min;
    Vector2 bb_max;

    Rectangle() : bb_min(Vector2()), bb_max(Vector2()) {}
    Rectangle(Vector2 bb_min, Vector2 bb_max) : bb_min(bb_min), bb_max(bb_max) {}

    Rectangle merge(const Rectangle& other) const {
        Vector2 new_bb_min = Vector2(std::min(bb_min.x, other.bb_min.x),
            std::min(bb_min.y, other.bb_min.y));
        Vector2 new_bb_max = Vector2(std::max(bb_max.x, other.bb_max.x),
            std::max(bb_max.y, other.bb_max.y));
        return Rectangle(new_bb_min, new_bb_max);
    }

    Vector2 center() const {
        return (bb_min + bb_max) / 2.0;
    }

    Vector2 size() const {
        return bb_max - bb_min;
    }

    bool intersects_ray(const Vector2& ray_origin, const Vector2& ray_direction) const {
        //Slab method
        //https://en.wikipedia.org/wiki/Slab_method
        Vector2 t_low = (bb_min - ray_origin) / ray_direction;
        Vector2 t_high = (bb_max - ray_origin) / ray_direction;
        Vector2 t_close = t_low.min(t_high);
        Vector2 t_far = t_low.max(t_high);

        real t_close_max = std::max(t_close.x, t_close.y);
        real t_far_min = std::min(t_far.x, t_far.y);
        return t_close_max <= t_far_min;
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
    Plane(const Vector3& normal, real dot_origin) : normal(normal), dot_origin(dot_origin) {}
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

//Vector3 tetrahedron_circumcenter(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3);
class Math {
public:
    static int wrap(int value, int min_val, int max_val) {
        int off_val = value - min_val;
        int range = max_val - min_val;
        return off_val < 0 ? off_val % range + range + min_val : off_val % range + min_val;
    }

    static real det(Vector2 a, Vector2 b) {
        return a.x * b.y - a.y * b.x;
    }

    static bool triangle_contains_point(const Vector2& p, const Vector2& p0, const Vector2& p1, const Vector2& p2) {
        //Barycentric coords
        real area = Math::det(p1 - p0, p2 - p0);
        real area_0 = Math::det(p1 - p, p2 - p);
        real area_1 = Math::det(p - p0, p2 - p0);
        real area_2 = Math::det(p1 - p0, p - p0);

        return signbit(area) == signbit(area_0) && signbit(area) == signbit(area_1) && signbit(area) == signbit(area_2);
    }

    static Vector2 triangle_circumcenter(const Vector2& p0, const Vector2& p1, const Vector2& p2) {
        return intersect_lines((p1 + p0) / 2.0, (p1 - p0).rot_CCW_90(), (p1 + p2) / 2.0, (p2 - p1).rot_CCW_90());
    }

    //@param p0 Point on line 0
    //@param r0 Ray pointing along line 0
    //@param p1 Point on line 1
    //@param r1 Ray pointing along line 1
    static Vector2 intersect_lines(const Vector2& p0, const Vector2& r0, const Vector2& p1, const Vector2& r1) {
        //Find [s, t] such that p0 + s * r0 == p1 + t * r1
        real det_r = det(r0, r1);
        if (det_r == 0)
            return p0;
            
        Vector2 dp = p1 - p0;
        real det_s = det(dp, r1);
        real s = det_s / det_r;
        return p0 + r0 * s;
    }

    static Vector3 tetrahedron_circumcenter(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) {
        //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter

        //From Matthias Muller
        //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68
        Vector3 b = p1 - p0;
        Vector3 c = p2 - p0;
        Vector3 d = p3 - p0;

        real det = 2.0 * (b.x * (c.y * d.z - c.z * d.y)
            - b.y * (c.x * d.z - c.z * d.x)
            + b.z * (c.x * d.y - c.y * d.x));

        if (det == 0.0) {
            return p0;
        }
        else {
            Vector3 v = c.cross(d) * b.dot(b) + d.cross(b) * c.dot(c) + b.cross(c) * d.dot(d);
            v /= det;
            return p0 + v;
        }
    }

};
    
}

#endif //CYCLOPS_MATH_H