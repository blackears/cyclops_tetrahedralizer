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

#include <iostream>

using namespace CyclopsTetra3D;

struct BVHTreeEdge2 {
    Vector2 p0;
    Vector2 p1;
    Vector2 centroid;
    Rectangle bounds;

    BVHTreeEdge2(const Vector2& p0, const Vector2& p1)
        : p0(p0), p1(p1), 
        bounds(Rectangle(p0.min(p1), p0.max(p1))),
        centroid((p0 + p1) / 2.0) {
    }

    bool intersect_ray(const Vector2& ray_origin, const Vector2& ray_direction, Vector2& out_hit_pos, Vector2& out_hit_normal) const {
        Vector2 line_dir = p1 - p0;
        real det_denom = Math::det(line_dir, ray_direction);
        if (det_denom == 0)
            return false;

        real det_s = Math::det(ray_origin - p0, ray_direction);

        real s = det_s / det_denom;
        if (s < 0 || s > 1)
            return false;

        out_hit_pos = p0 + line_dir * s;
        out_hit_normal = (ray_origin - out_hit_pos).normalized();

        return true;
    }
};

struct BVHTreeNode2 {
    Rectangle bounds;

    unsigned int child_left_idx;
    //Right child is child_left_idx + 1

    unsigned int first_edge_offset;
    unsigned int num_edges;

    bool is_leaf() const {
        return num_edges > 0;
    }

    void update_bounds(const std::vector<BVHTreeEdge2>& edges, const std::vector<BVHTreeNode2>& nodes) {
        if (num_edges > 0) {
            //leaf node with edges as children
            bounds = edges[first_edge_offset].bounds;
            for (unsigned int i = 1; i < num_edges; i++) {
                bounds = bounds.merge(edges[first_edge_offset + i].bounds);
            }
        }
        else {
            //intermediate node with two nodes as children
            bounds = nodes[child_left_idx].bounds.merge(nodes[child_left_idx + 1].bounds);
        }
    }

    bool ray_cast(const Vector2& ray_origin, const Vector2& ray_direction,
        const std::vector<BVHTreeEdge2>& edges, const std::vector<BVHTreeNode2>& nodes,
        Vector2& out_hit_pos, Vector2& out_hit_normal, int& out_index) const {

        std::cout << "ray_cast: bounds " << bounds << std::endl;

        if (!bounds.intersects_ray(ray_origin, ray_direction)) {
            std::cout << "  missed bounds" << std::endl;
            return false;
        }

        if (is_leaf()) {
            for (unsigned int i = 0; i < num_edges; i++) {
                if (edges[first_edge_offset + i].intersect_ray(ray_origin, ray_direction, out_hit_pos, out_hit_normal)) {
                    std::cout << "  isect edge" << std::endl;
                    out_index = first_edge_offset + i;
                    return true;
                }
            }
            return false;
        }
        else {
            return nodes[child_left_idx].ray_cast(ray_origin, ray_direction, edges, nodes, out_hit_pos, out_hit_normal, out_index) ||
                nodes[child_left_idx + 1].ray_cast(ray_origin, ray_direction, edges, nodes, out_hit_pos, out_hit_normal, out_index);
        }
    }

};

class BVHTree2 {
    std::vector<BVHTreeEdge2> edges;
    std::vector<BVHTreeNode2> nodes;

    void build_nodes_recursive(int root_node_idx, int& nodes_used);

    static const Vector2 face_check_dirs[];

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    void build_from_edges(const std::vector<Vector2>& points, const std::vector<int>& indices);

    bool ray_cast(const Vector2& ray_origin, const Vector2& ray_direction, Vector2& out_hit_pos, Vector2& out_hit_normal, int& out_index) const;

    bool is_inside(const Vector2& p, real epsilon = 0.0) const;

};

#endif //CYCLOPS_BVH_TREE2_H