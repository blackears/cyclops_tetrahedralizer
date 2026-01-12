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

/*
 * Substantial portions of this are based on Matthias Muller's tetrahedralizer adon for Blender.
 * https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py
 */

#ifndef CYCLOPS_TESSELLATE_H
#define CYCLOPS_TESSELLATE_H

#include  <vector>

namespace CyclopsTessellate3D {

//typedef float real_t;

template<typename T = float>
struct Vector3 {
    T x;
    T y;
    T z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    T magnitude() const { return sqrt(x * x + y * y + z * z); }
    Vector3<T> cross(const Vector3<T> rhs) const { return Vector3<T>(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x); }
    T dot(const Vector3<T> rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }

    Vector3<T> min(const Vector3<T> rhs) const { return Vector3<T>(min(x, rhs.x), min(y, rhs.y), min(z, rhs.z)); }
    Vector3<T> max(const Vector3<T> rhs) const { return Vector3<T>(max(x, rhs.x), max(y, rhs.y), max(z, rhs.z)); }

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

    Vector3<T>& operator*=(T rhs) {
        this->x *= rhs;
        this->y *= rhs;
        return *this;
    }

    friend Vector3<T> operator*(Vector3<T> lhs, T rhs) {
        lhs *= rhs;
        return lhs;
    }

    Vector3<T>& operator/=(T rhs) {
        this->x /= rhs;
        this->y /= rhs;
        return *this;
    }

    friend Vector3<T> operator/(Vector3<T> lhs, T rhs) {
        lhs /= rhs;
        return lhs;
    }
};

// struct Triangle {
//     int v0_idx;
//     int v1_idx;
//     int v2_idx;
// };

template<typename T = float>
struct Tetrahedron {
    int v0_idx;
    int v1_idx;
    int v2_idx;
    int v3_idx;

    Vector3<T> calc_circum_center(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const;

    T quality(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const;
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

};


template<typename T = float>
struct BVHTreeTriangle {
    Vector3<T> p0;
    Vector3<T> p1;
    Vector3<T> p2; 
    Vector3<T> centroid;
    BoundingBox<T> bounds;

    BVHTreeTriangle(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2)
        : p0(p0), p1(p1), p2(p2),
        bounds(BoundingBox<T>(p0.min(p1).min(p2), p0.max(p1).max(p2))),
        centroid((p0 + p1 + p2) / 3.0) {}
};

template<typename T = float>
struct BVHTreeNode {
    BoundingBox<T> bounds;

    unsigned int child_left_idx;
    //Right child is child_left_idx + 1

    unsigned int first_primitive_offset;
    unsigned int num_primitives;

    void update_bounds(const std::vector<BVHTreeTriangle<T>>& primitives, const std::vector<BVHTreeNode<T>>& nodes) {
        if (num_primitives > 0) {
            bounds = primitives[first_primitive_offset].bounds;
            for (unsigned int i = 1; i < num_primitives; i++) {
                bounds = bounds.merge(primitives[first_primitive_offset + i].bounds);
            }
        } else {
            bounds = nodes[child_left_idx].bounds.merge(nodes[child_left_idx + 1].bounds);
        }
    }
};

//https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
template<typename T = float>
class BVHTree {
    std::vector<BVHTreeTriangle<T>> triangles;
    std::vector<BVHTreeNode<T>> nodes;

    void build_nodes_recursive(int root_node_idx, int& nodes_used);

public:
    void build_from_triangles(const std::vector<Vector3<T>>& points, const std::vector<int>& indices);

    void ray_cast(Vector3<T> origin, Vector3<T> direction, T distance, Vector3<T> &out_hit_pos, Vector3<T> &out_hit_normal, int &out_index, T &out_hit_distance) const;

    bool is_inside(Vector3<T> p, T dist_min = 0.0) const;

};



template<typename T = float>
class CyclopsTetrahedralizer {
    Vector3<T>* point_list;

private:
    void create_tetrahedron_ids(const std::vector<Vector3<T>>& points, BVHTree<T>& bvh_tree, float quality_threshold);

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    //@param resolution spacing for extra interior points
    //@param quality_threshold minimum quality for tetrahedrons.  0.0 to 1.0, with 1.0 being equilateral tetrahedra.
    void create_tetrahedrons(const std::vector<Vector3<T>>& points, const std::vector<int>& indices, float resolution = .1, float quality_threshold = 0.001);


    //void tessellate_tetrahedra(float* mesh_points);

};

}

#endif //CYCLOPS_TESSELLATE_H