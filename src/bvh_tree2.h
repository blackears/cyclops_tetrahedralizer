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

#ifndef CYCLOPS_BVH_TREE2_H
#define CYCLOPS_BVH_TREE2_H

#include "math.h"

using namespace CyclopsTetra3D;

struct BVHTreeTriangle2 {
    Vector2 p0;
    Vector2 p1;
    Vector2 p2;
    Vector2 centroid;
    Rectangle bounds;

    BVHTreeTriangle2(const Vector2& p0, const Vector2& p1, const Vector2& p2)
        : p0(p0), p1(p1), p2(p2),
        bounds(Rectangle(p0.min(p1).min(p2), p0.max(p1).max(p2))),
        centroid((p0 + p1 + p2) / 3.0) {
    }

    bool intersect_ray(const Vector2& ray_origin, const Vector2& ray_direction, Vector2& out_hit_pos, Vector2& out_hit_normal) const {
        real dir_sq = ray_direction.dot(ray_direction);
        real p0_hit = (p0 - ray_origin).dot(ray_direction) / dir_sq;
        real p1_hit = (p1 - ray_origin).dot(ray_direction) / dir_sq;
        real p2_hit = (p2 - ray_origin).dot(ray_direction) / dir_sq;

        if ((p0_hit <= 1.0 && p1_hit >= 1.0) || (p0_hit >= 1.0 && p1_hit <= 1.0))
            return true;
        if ((p0_hit <= 1.0 && p2_hit >= 1.0) || (p0_hit >= 1.0 && p2_hit <= 1.0))
            return true;
        if ((p1_hit <= 1.0 && p2_hit >= 1.0) || (p1_hit >= 1.0 && p2_hit <= 1.0))
            return true;

        return false;
    }
};

struct BVHTreeNode2 {
    Rectangle bounds;

    unsigned int child_left_idx;
    //Right child is child_left_idx + 1

    unsigned int first_triangle_offset;
    unsigned int num_triangles;

    bool is_leaf() const {
        return num_triangles > 0;
    }

    void update_bounds(const std::vector<BVHTreeTriangle2>& triangles, const std::vector<BVHTreeNode2>& nodes) {
        if (num_triangles > 0) {
            bounds = triangles[first_triangle_offset].bounds;
            for (unsigned int i = 1; i < num_triangles; i++) {
                bounds = bounds.merge(triangles[first_triangle_offset + i].bounds);
            }
        }
        else {
            bounds = nodes[child_left_idx].bounds.merge(nodes[child_left_idx + 1].bounds);
        }
    }

    bool ray_cast(const Vector2& ray_origin, const Vector2& ray_direction,
        const std::vector<BVHTreeTriangle2>& triangles, const std::vector<BVHTreeNode2>& nodes,
        Vector2& out_hit_pos, Vector2& out_hit_normal, int& out_index) const {

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
        }
        else {
            return nodes[child_left_idx].ray_cast(ray_origin, ray_direction, triangles, nodes, out_hit_pos, out_hit_normal, out_index) ||
                nodes[child_left_idx + 1].ray_cast(ray_origin, ray_direction, triangles, nodes, out_hit_pos, out_hit_normal, out_index);
        }
    }

};

class BVHTree2 {
    std::vector<BVHTreeTriangle2> triangles;
    std::vector<BVHTreeNode2> nodes;

    void build_nodes_recursive(int root_node_idx, int& nodes_used);

    static const Vector2 face_check_dirs[];

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    void build_from_triangles(const std::vector<Vector2>& points, const std::vector<int>& indices);

    bool ray_cast(const Vector2& ray_origin, const Vector2& ray_direction, Vector2& out_hit_pos, Vector2& out_hit_normal, int& out_index) const;

    bool is_inside(const Vector2& p, real dist_min = 0.0) const;

};

#endif //CYCLOPS_BVH_TREE2_H