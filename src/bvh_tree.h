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
  * BVH Tree implementation based on jbikker's blog series:
  * https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
  */

#ifndef CYCLOPS_BVH_TREE_H
#define CYCLOPS_BVH_TREE_H

#include "math.h"

using namespace CyclopsTetra3D;
 
struct BVHTreeTriangle {
    Vector3 p0;
    Vector3 p1;
    Vector3 p2; 
    Vector3 centroid;
    BoundingBox bounds;

    BVHTreeTriangle(const Vector3& p0, const Vector3& p1, const Vector3& p2)
        : p0(p0), p1(p1), p2(p2),
        bounds(BoundingBox(p0.min(p1).min(p2), p0.max(p1).max(p2))),
        centroid((p0 + p1 + p2) / 3.0) {}
    
    bool intersect_ray(const Vector3& ray_origin, const Vector3& ray_direction, Vector3 &out_hit_pos, Vector3 &out_hit_normal) const {
        Vector3 edge1 = p1 - p0;
        Vector3 edge2 = p2 - p0;
        Vector3 h = ray_direction.cross(edge2);
        real a = edge1.dot(h);
        if (a > -1e-5 && a < 1e-5) {
            return false; // This ray is parallel to this triangle.
        }
        real f = 1.0 / a;
        Vector3 s = ray_origin - p0;
        real u = f * s.dot(h);
        if (u < 0.0 || u > 1.0) {
            return false;
        }
        Vector3 q = s.cross(edge1);
        real v = f * ray_direction.dot(q);
        if (v < 0.0 || u + v > 1.0) {
            return false;
        }
        // At this stage we can compute t to find out where the intersection point is on the line.
        real t = f * edge2.dot(q);
        if (t > 1e-5) { // ray intersection
            out_hit_pos = ray_origin + ray_direction * t;
            out_hit_normal = edge1.cross(edge2).normalized();
            return true;
        }
        
        // This means that there is a line intersection but not a ray intersection.
        return false;
    }
};

struct BVHTreeNode {
    BoundingBox bounds;

    unsigned int child_left_idx;
    //Right child is child_left_idx + 1

    unsigned int first_triangle_offset;
    unsigned int num_triangles;

    bool is_leaf() const {
        return num_triangles > 0;
    }

    void update_bounds(const std::vector<BVHTreeTriangle>& triangles, const std::vector<BVHTreeNode>& nodes) {
        if (num_triangles > 0) {
            bounds = triangles[first_triangle_offset].bounds;
            for (unsigned int i = 1; i < num_triangles; i++) {
                bounds = bounds.merge(triangles[first_triangle_offset + i].bounds);
            }
        } else {
            bounds = nodes[child_left_idx].bounds.merge(nodes[child_left_idx + 1].bounds);
        }
    }

    bool ray_cast(const Vector3& ray_origin, const Vector3& ray_direction, 
        const std::vector<BVHTreeTriangle>& triangles, const std::vector<BVHTreeNode>& nodes,
        Vector3 &out_hit_pos, Vector3 &out_hit_normal, int &out_index) const {

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

class BVHTree {
    std::vector<BVHTreeTriangle> triangles;
    std::vector<BVHTreeNode> nodes;

    void build_nodes_recursive(int root_node_idx, int& nodes_used);

    static const Vector3 face_check_dirs[];

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    void build_from_triangles(const std::vector<Vector3>& points, const std::vector<int>& indices);

    bool ray_cast(const Vector3& ray_origin, const Vector3& ray_direction, Vector3 &out_hit_pos, Vector3 &out_hit_normal, int &out_index) const;

    bool is_inside(const Vector3& p, real dist_min = 0.0) const;

};

#endif //CYCLOPS_BVH_TREE_H