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

#include  <vector>

#include "math.h"

namespace CyclopsTessellate3D {

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
    
    bool intersect_ray(const Vector3<T>& ray_origin, const Vector3<T>& ray_direction, Vector3<T> &out_hit_pos, Vector3<T> &out_hit_normal) const {
        Vector3<T> edge1 = p1 - p0;
        Vector3<T> edge2 = p2 - p0;
        Vector3<T> h = ray_direction.cross(edge2);
        T a = edge1.dot(h);
        if (a > -1e-5 && a < 1e-5) {
            return false; // This ray is parallel to this triangle.
        }
        T f = 1.0 / a;
        Vector3<T> s = ray_origin - p0;
        T u = f * s.dot(h);
        if (u < 0.0 || u > 1.0) {
            return false;
        }
        Vector3<T> q = s.cross(edge1);
        T v = f * ray_direction.dot(q);
        if (v < 0.0 || u + v > 1.0) {
            return false;
        }
        // At this stage we can compute t to find out where the intersection point is on the line.
        T t = f * edge2.dot(q);
        if (t > 1e-5) { // ray intersection
            out_hit_pos = ray_origin + ray_direction * t;
            out_hit_normal = edge1.cross(edge2).normalize();
            return true;
        }
        
        // This means that there is a line intersection but not a ray intersection.
        return false;
    }
};

template<typename T = float>
struct BVHTreeNode {
    BoundingBox<T> bounds;

    unsigned int child_left_idx;
    //Right child is child_left_idx + 1

    unsigned int first_triangle_offset;
    unsigned int num_triangles;

    bool is_leaf() const {
        return num_triangles > 0;
    }

    void update_bounds(const std::vector<BVHTreeTriangle<T>>& triangles, const std::vector<BVHTreeNode<T>>& nodes) {
        if (num_triangles > 0) {
            bounds = triangles[first_triangle_offset].bounds;
            for (unsigned int i = 1; i < num_triangles; i++) {
                bounds = bounds.merge(triangles[first_triangle_offset + i].bounds);
            }
        } else {
            bounds = nodes[child_left_idx].bounds.merge(nodes[child_left_idx + 1].bounds);
        }
    }

    bool ray_cast(Vector3<T> ray_origin, Vector3<T> ray_direction, 
        const std::vector<BVHTreeTriangle<T>>& triangles, const std::vector<BVHTreeNode<T>>& nodes,
        Vector3<T> &out_hit_pos, Vector3<T> &out_hit_normal, int &out_index) const {

        if (!bounds.intersects_ray(ray_origin, ray_direction))
            return false;
        if (is_leaf()) {
            for (unsigned int i = 0; i < num_triangles; i++) {
                if (triangles[first_triangle_offset + i].intersect_ray(ray_origin, ray_direction, out_hit_pos, out_hit_normal)) {
                    out_index = first_triangle_offset + i;
                    return true;
                }
            }
            return false;
        } else {
            return nodes[child_left_idx].ray_cast(ray_origin, ray_direction, triangles, nodes, out_hit_pos, out_hit_normal, out_index) ||
                nodes[child_left_idx + 1].ray_cast(ray_origin, ray_direction, triangles, nodes, out_hit_pos, out_hit_normal, out_index);
        }
    }

    // bool intersects_ray(Vector3<T> ray_origin, Vector3<T> ray_direction, const std::vector<BVHTreeTriangle<T>>& primitives, const std::vector<BVHTreeNode<T>>& nodes) {
    //     if (bounds.intersects_ray(ray_origin, ray_direction))
    //         return true;
    //     if (is_leaf()) {
    //         for (unsigned int i = 0; i < num_triangles; i++) {
    //             if (primitives[first_triangle_offset + i].bounds.intersects_ray(ray_origin, ray_direction))
    //                 return true;
    //         }
    //         return false;
    //     } else {
    //         return nodes[child_left_idx].intersects_ray(ray_origin, ray_direction, primitives, nodes) ||
    //             nodes[child_left_idx + 1].intersects_ray(ray_origin, ray_direction, primitives, nodes);
    //     }
    // }

};

//https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
template<typename T = float>
class BVHTree {
    std::vector<BVHTreeTriangle<T>> triangles;
    std::vector<BVHTreeNode<T>> nodes;

    void build_nodes_recursive(int root_node_idx, int& nodes_used);

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    void build_from_triangles(const std::vector<Vector3<T>>& points, const std::vector<int>& indices);

    bool ray_cast(Vector3<T> ray_origin, Vector3<T> ray_direction, Vector3<T> &out_hit_pos, Vector3<T> &out_hit_normal, int &out_index) const;

    bool is_inside(Vector3<T> p, T dist_min = 0.0) const;

};

template<typename T = float>
struct Plane {
    Vector3<T> normal;
    //distance along normal from origin to plane
    T dot_origin;

    Plane() : normal(Vector3<T>()), dot_origin(0) {}
    Plane(const Vector3<T>& normal, T origin) : normal(normal), dot_origin(origin.dot(normal)) {}
    Plane(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2) : normal((p1 - p0).cross(p2 - p0)), dot_origin(p0.dot(normal)) {}
    
    T distance_to_plane(const Vector3<T>& p) const {
        return normal.dot(p) - dot_origin;
    }

    bool intersect_ray(const Vector3<T>& ray_origin, const Vector3<T>& ray_direction, Vector3<T>& out_intersection) const {
        T denom = normal.dot(ray_direction);
        if (denom == 0.0) {
            return false;
        }
        T numer = dot_origin - normal.dot(ray_origin);

        T s = numer / denom;
        out_intersection = ray_origin + ray_direction * s;
        return true;
    }
};


//Face winding - face normals points outward
template<typename T = float>
struct Tetrahedron {
    //Vertex ordering per face
    static constexpr int face_vert_indices[4][3] = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {0, 2, 3}};
    static constexpr int face_missing_vert_index[4] = {3, 2, 0, 1};

    int vert_indices[4];
    int neighbors[4];

    Vector3<T> circumcenter;
    Vector3<T> center;

    Plane<T> face_planes[4];
    //bool boundary_face[4];

    bool valid;

    void create_from_points(int v0_idx, int v1_idx, int v2_idx, int v3_idx, const std::vector<Vector3<T>>& points);

    bool has_neighbor(int neighbor_idx) const {
        for (int i = 0; i < 4; i++) {
            if (neighbors[i] == neighbor_idx) {
                return true;
            }
        }
        return false;
    }
    Vector3<T> calc_circumcenter(const std::vector<Vector3<T>>& points) const;
    bool point_in_circumsphere(const Vector3<T>& p, const std::vector<Vector3<T>>& points) const {
        return (circumcenter - p).magnitude_squared() < (circumcenter - points[v0_idx]).magnitude_squared();
    }
    bool contains_point(const Vector3<T>& p, const std::vector<Vector3<T>>& points) const;
    int find_adjacent_tetrahedron(const Vector3<T>& dir, const std::vector<Vector3<T>>& points) const;
    T quality(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& p3) const;

    //vertex indices must wind face outward
    int find_face(int v0_idx, int v1_idx, int v2_idx) const {
        for (int i = 0; i < 4; i++) {
            if ((vert_indices[face_vert_indices[i][0]] == v0_idx &&
                 vert_indices[face_vert_indices[i][1]] == v1_idx &&
                 vert_indices[face_vert_indices[i][2]] == v2_idx) ||
                (vert_indices[face_vert_indices[i][0]] == v1_idx &&
                 vert_indices[face_vert_indices[i][1]] == v2_idx &&
                 vert_indices[face_vert_indices[i][2]] == v0_idx) ||
                (vert_indices[face_vert_indices[i][0]] == v2_idx &&
                 vert_indices[face_vert_indices[i][1]] == v0_idx &&
                 vert_indices[face_vert_indices[i][2]] == v1_idx)) {
                return i;
            }
        }
        return -1;

    }
};

template<typename T = float>
class CyclopsTetrahedralizer {
    Vector3<T>* point_list;

private:
//    void create_tetrahedrons_internal(const std::vector<Vector3<T>>& points, BVHTree<T>& bvh_tree, float quality_threshold);
    void create_tetrahedrons_iter(std::vector<Tetrahedron<T>>& tetrahedrons, const std::vector<Vector3<T>>& points);

//    Vector3<T> calc_circum_center(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const;

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