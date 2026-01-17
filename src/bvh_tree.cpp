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

#include "bvh_tree.h"

#include <algorithm>

using namespace CyclopsTetra3D;


void BVHTree::build_nodes_recursive(int root_node_idx, int& nodes_used) {
    BVHTreeNode& root_node = nodes[root_node_idx];

    Vector3 bb_size = root_node.bounds.size();
    int axis;
    real midpoint;
    if (bb_size.x >= bb_size.y && bb_size.x >= bb_size.z) {
        axis = 0;
        midpoint = bb_size.x / 2.0 + root_node.bounds.bb_min.x;
    }
    else if (bb_size.y >= bb_size.z) {
        axis = 1;
        midpoint = bb_size.y / 2.0 + root_node.bounds.bb_min.y;
    }
    else {
        axis = 2;
        midpoint = bb_size.z / 2.0 + root_node.bounds.bb_min.z;
    }

    //Partition primitives around midpoint - sort subpartions
    int i = root_node.first_triangle_offset;
    int j = i + root_node.num_triangles - 1;
    while (i <= j) {
        BVHTreeTriangle& tri = triangles[i];
        if (tri.centroid[axis] < midpoint) {
            i++;
        } else {
            std::swap(triangles[i], triangles[j]);
            j--;
        }
    }
    
    int left_count = i - root_node.first_triangle_offset;
    if (left_count == 0 || left_count == root_node.num_triangles) {
        return; //Cannot split further
    }

    //Next two nodes in node list with be children of this node
    nodes[nodes_used].first_triangle_offset = root_node.first_triangle_offset;
    nodes[nodes_used].num_triangles = left_count;
    nodes[nodes_used].update_bounds(triangles, nodes);

    nodes[nodes_used + 1].first_triangle_offset = i;
    nodes[nodes_used + 1].num_triangles = root_node.num_triangles - left_count;
    nodes[nodes_used + 1].update_bounds(triangles, nodes);
    
    root_node.child_left_idx = nodes_used;
    root_node.num_triangles = 0;
    nodes_used += 2;

    //Recurse
    build_nodes_recursive(root_node.child_left_idx, nodes_used);
    build_nodes_recursive(root_node.child_left_idx + 1, nodes_used);
}

void BVHTree::build_from_triangles(const std::vector<Vector3>& points, const std::vector<int>& indices) {
    nodes.clear();
    triangles.clear();

    nodes.resize(indices.size() / 3 * 2 - 1);
    triangles.reserve(indices.size() / 3);

    for (int i = 0; i < indices.size() / 3; i++) {
        int idx0 = indices[i * 3 + 0];
        int idx1 = indices[i * 3 + 1];
        int idx2 = indices[i * 3 + 2];

        Vector3 v0 = points[idx0];
        Vector3 v1 = points[idx1];
        Vector3 v2 = points[idx2];

        BVHTreeTriangle tri = BVHTreeTriangle(v0, v1, v2);
        triangles.push_back(tri);
    }

    nodes[0].child_left_idx = 0;
    //nodes[0].right_child_idx = 0;
    nodes[0].first_triangle_offset = 0;
    nodes[0].num_triangles = triangles.size();
    nodes[0].update_bounds(triangles, nodes);

    int nodes_used = 1;
    build_nodes_recursive(0, nodes_used);
}

bool BVHTree::is_inside(const Vector3& p, real dist_min) const {
    //Check multiple directions to minimize chance of hitting an edge

    int num_inside = 0;

    for (int i = 0; i < 6; i++) {
        Vector3 hit_pos;
        Vector3 hit_normal;
        int hit_index;

        nodes[0].ray_cast(p, face_check_dirs[i], triangles, nodes, hit_pos, hit_normal, hit_index);
        real hit_distance = (hit_pos - p).magnitude();
        //Check hit is valid
        if (hit_index >= 0) {
            if (hit_normal.dot(face_check_dirs[i]) > 0) {
                num_inside++;
            }

            if (dist_min > 0.0 && hit_distance < dist_min) {
                //If hit outside face, we are outside
                return false;
            }
        }
    }

    return num_inside >= 3;
}

bool BVHTree::ray_cast(const Vector3& ray_origin, const Vector3& ray_direction, Vector3 &hit_pos, Vector3 &hit_normal, int &out_index) const {
    return nodes[0].ray_cast(ray_origin, ray_direction, triangles, nodes, hit_pos, hit_normal, out_index);
}

