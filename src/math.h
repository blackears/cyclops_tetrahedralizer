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

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "predicates/predicates.h"

namespace CyclopsTetra3D {

typedef double real;
//typedef float real;


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

    Vector2 rot_CCW_90() const { return Vector2(-y, x); }
    Vector2 rot_CW_90() const { return Vector2(y, -x); }
    Vector2 reverse() const { return Vector2(y, x); }

    std::string to_string() const { return "(" + std::to_string(x) + ", " + std::to_string(y) + ")"; }

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


    Vector2& operator*=(const Vector2& rhs) {
        this->x *= rhs.x;
        this->y *= rhs.y;
        return *this;
    }

    friend Vector2 operator*(Vector2 lhs, const Vector2& rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector2& operator/=(const Vector2& rhs) {
        this->x /= rhs.x;
        this->y /= rhs.y;
        return *this;
    }

    friend Vector2 operator/(Vector2 lhs, const Vector2& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const Vector2& lhs, const Vector2& rhs) {
        return lhs.x == rhs.x && lhs.y == rhs.y;
    }

    friend bool operator<(const Vector2& lhs, const Vector2& rhs) {
        if (lhs.x != rhs.x)
            return lhs.x < rhs.x;
        return lhs.y < rhs.y;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector2& obj) {
        os << "(" << obj.x << ", " << obj.y << ")";
        return os;
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


    Vector3& operator*=(const Vector3& rhs) {
        this->x *= rhs.x;
        this->y *= rhs.y;
        this->z *= rhs.z;
        return *this;
    }

    friend Vector3 operator*(Vector3 lhs, const Vector3& rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3& operator/=(const Vector3& rhs) {
        this->x /= rhs.x;
        this->y /= rhs.y;
        this->z /= rhs.z;
        return *this;
    }

    friend Vector3 operator/(Vector3 lhs, const Vector3& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const Vector3& lhs, const Vector3& rhs) {
        return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
    }

    friend bool operator<(const Vector3& lhs, const Vector3& rhs) {
        if (lhs.x != rhs.x)
            return lhs.x < rhs.x;
        if (lhs.y != rhs.y)
            return lhs.y < rhs.y;
        return lhs.z < rhs.z;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3& obj) {
        os << "(" << obj.x << ", " << obj.y << ", " << obj.z << ")";
        return os;
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

    friend std::ostream& operator<<(std::ostream& os, const Rectangle& obj) {
        os << "(" << obj.bb_min << ", " << obj.bb_max << ")";
        return os;
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

    bool intersects_ray(const Vector3& ray_origin, const Vector3& ray_direction, real epsilon = 1e-4) const {
        Vector3 eps_vec(epsilon, epsilon, epsilon);
        //Slab method
        //https://en.wikipedia.org/wiki/Slab_method
        Vector3 t_low = (bb_min - eps_vec - ray_origin) / ray_direction;
        Vector3 t_high = (bb_max + eps_vec - ray_origin) / ray_direction;
        Vector3 t_close = t_low.min(t_high);
        Vector3 t_far = t_low.max(t_high);

        real t_close_max = std::max(std::max(t_close.x, t_close.y), t_close.z);
        real t_far_min = std::min(std::min(t_far.x, t_far.y), t_far.z);
        return t_close_max <= t_far_min;
    }

    friend std::ostream& operator<<(std::ostream& os, const BoundingBox& obj) {
        os << "(" << obj.bb_min << ", " << obj.bb_max << ")";
        return os;
    }
};


struct Plane {
    Vector3 normal;
    //distance along normal from origin to plane
    real dot_origin;

    Plane() : normal(Vector3()), dot_origin(0) {}
    Plane(const Vector3& normal, real dot_origin) : normal(normal), dot_origin(dot_origin) {}

    static Plane create(const Vector3& normal, const Vector3& p) {
        return Plane(normal, p.dot(normal));
    }

    static Plane create(const Vector3& p0, const Vector3& p1, const Vector3& p2) {
        Vector3 normal = ((p1 - p0).cross(p2 - p0)).normalized();
        real dot = p0.dot(normal);
        return Plane(normal, dot);
    }

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

    friend std::ostream& operator<<(std::ostream& os, const Plane& obj) {
        os << "(" << obj.normal << ", " << obj.dot_origin << ")";
        return os;
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

    // Determinant of a 2x2 matrix where [a b] are columns
    static real det(Vector2 a, Vector2 b) {
        return a.x * b.y - b.x * a.y;
    }

    // Determinant of a 3x3 matrix where [a b c] are columns
    static real det(Vector3 a, Vector3 b, Vector3 c) {
        return a.x * (b.y * c.z - c.y * b.z)
            + b.x * (c.y * a.z - a.y * c.z)
            + c.x * (a.y * b.z - b.y * a.z);
    }

    static real dist_to_segment_squared(const Vector3& p, const Vector3& p0, const Vector3& p1) {
        Vector3 a = p - p0;
        Vector3 b = p1 - p0;

        //Scalar for vector b that is projection of p onto segment [p1 - p0]
        real s = std::clamp(a.dot(b) / b.dot(b), real(0.0), real(1.0));

        return (b * s - p).magnitude_squared();
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

    static bool tetrahedron_contains_point(const Vector3& p, const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) {
        //Barycentric coords - these subregions should have the same signed area as the whole
        real area = Math::det(p1 - p0, p2 - p0, p3 - p0);
        real area_0 = Math::det(p1 - p, p2 - p, p3 - p);
        real area_1 = Math::det(p - p0, p2 - p0, p3 - p0);
        real area_2 = Math::det(p1 - p0, p - p0, p3 - p0);
        real area_3 = Math::det(p1 - p0, p2 - p0, p - p0);

        return signbit(area) == signbit(area_0) && signbit(area) == signbit(area_1) && signbit(area) == signbit(area_2) && signbit(area) == signbit(area_3);
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

    //@return a positive value if the points pa, pb, and pc occur in counterclockwise order; 
    // a negative value if they occur in clockwise order; and zero if they are collinear.  
    // The result is also a rough approximation of twice the signed area of the triangle 
    // defined by the three points.
    //
    // This uses exact arithmetic to ensure a correct answer and avoid round off errors for 
    // determinants close to zero.
    static real test_orient_2d(const Vector2& pa, const Vector2& pb, const Vector2& pc) {
        double a[2] = { pa.x, pa.y };
        double b[2] = { pb.x, pb.y };
        double c[2] = { pc.x, pc.y };

        return orient2d(a, b, c);
    }

    //@return a positive value if the point pd lies below the plane passing through 
    // pa, pb, and pc; "below" is defined so that pa, pb, and pc appear in 
    // counterclockwise order when viewed from above the plane.  Returns a negative 
    // value if pd lies above the plane.  Returns zero if the points are coplanar.  
    // The result is also a rough approximation of six times the signed volume of the 
    // tetrahedron defined by the four points.
    //
    // This uses exact arithmetic to ensure a correct answer and avoid round off errors for 
    // determinants close to zero.
    static real test_orient_3d(const Vector3& pa, const Vector3& pb, const Vector3& pc, const Vector3& pd) {
        double a[3] = { pa.x, pa.y, pa.z };
        double b[3] = { pb.x, pb.y, pb.z };
        double c[3] = { pc.x, pc.y, pc.z };
        double d[3] = { pd.x, pd.y, pd.z };

        return orient3d(a, b, c, d);
    }

    //@return a positive value if the point pd lies inside the circle passing 
    // through pa, pb, and pc; a negative value if it lies outside; and zero 
    // if the four points are cocircular.  The points pa, pb, and pc must be 
    // in counterclockwise order, or the sign of the result will be reversed.
    //
    // This uses exact arithmetic to ensure a correct answer and avoid round off errors for 
    // determinants close to zero.
    static real test_in_circle(const Vector2& pa, const Vector2& pb, const Vector2& pc, const Vector2& pd) {
        double a[2] = { pa.x, pa.y };
        double b[2] = { pb.x, pb.y };
        double c[2] = { pc.x, pc.y };
        double d[2] = { pd.x, pd.y };

        return incircle(a, b, c, d);
    }

    //@return a positive value if the point pe lies inside the sphere passing through 
    // pa, pb, pc, and pd; a negative value if it lies outside; and zero if the five 
    // points are cospherical.  The points pa, pb, pc, and pd must be ordered so that 
    // they have a positive orientation (as defined by orient3d()), or the sign of the 
    // result will be reversed.
    //
    // This uses exact arithmetic to ensure a correct answer and avoid round off errors for 
    // determinants close to zero.
    static real test_in_sphere(const Vector3& pa, const Vector3& pb, const Vector3& pc, const Vector3& pd, const Vector3& pe) {
        double a[3] = { pa.x, pa.y, pa.z };
        double b[3] = { pb.x, pb.y, pb.z };
        double c[3] = { pc.x, pc.y, pc.z };
        double d[3] = { pd.x, pd.y, pd.z };
        double e[3] = { pe.x, pe.y, pe.z };

        return insphere(a, b, c, d, e);
    }

};
    
}

#endif //CYCLOPS_MATH_H