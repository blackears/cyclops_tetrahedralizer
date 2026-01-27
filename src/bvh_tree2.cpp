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

#include "bvh_tree2.h"

#include <algorithm>

using namespace CyclopsTetra3D;


const Vector2 BVHTree2::face_check_dirs[] = {
    Vector2(1, 0),
    Vector2(-1, 0),
    Vector2(0, 1),
    Vector2(0, -1),
};

void BVHTree2::build_nodes_recursive(int root_node_idx, int& nodes_used) {
    BVHTreeNode2& root_node = nodes[root_node_idx];

    Vector2 bb_size = root_node.bounds.size();
    int axis;
    real midpoint;
    if (bb_size.x >= bb_size.y) {
        axis = 0;
        midpoint = bb_size.x / 2.0 + root_node.bounds.bb_min.x;
    }
    else {
        axis = 1;
        midpoint = bb_size.y / 2.0 + root_node.bounds.bb_min.y;
    }

    //Partition primitives around midpoint - sort subpartions
    int i = root_node.first_edge_offset;
    int j = i + root_node.num_edges - 1;
    while (i <= j) {
        BVHTreeEdge2& tri = edges[i];
        if (tri.centroid[axis] < midpoint) {
            i++;
        }
        else {
            std::swap(edges[i], edges[j]);
            j--;
        }
    }

    int left_count = i - root_node.first_edge_offset;
    if (left_count == 0 || left_count == root_node.num_edges) {
        return; //Cannot split further
    }

    //Next two nodes in node list with be children of this node
    nodes[nodes_used].first_edge_offset = root_node.first_edge_offset;
    nodes[nodes_used].num_edges = left_count;
    nodes[nodes_used].update_bounds(edges, nodes);

    nodes[nodes_used + 1].first_edge_offset = i;
    nodes[nodes_used + 1].num_edges = root_node.num_edges - left_count;
    nodes[nodes_used + 1].update_bounds(edges, nodes);

    root_node.child_left_idx = nodes_used;
    root_node.num_edges = 0;
    nodes_used += 2;

    //Recurse
    build_nodes_recursive(root_node.child_left_idx, nodes_used);
    build_nodes_recursive(root_node.child_left_idx + 1, nodes_used);
}

void BVHTree2::build_from_edges(const std::vector<Vector2>& points, const std::vector<int>& indices) {
    nodes.clear();
    edges.clear();

    nodes.resize(indices.size() - 1);
    edges.reserve(indices.size() / 2);

    for (int i = 0; i < indices.size() / 2; i++) {
        int idx0 = indices[i * 2 + 0];
        int idx1 = indices[i * 2 + 1];

        Vector2 v0 = points[idx0];
        Vector2 v1 = points[idx1];

        edges.push_back(BVHTreeEdge2(v0, v1));
    }

    nodes[0].child_left_idx = 0;
    //nodes[0].right_child_idx = 0;
    nodes[0].first_edge_offset = 0;
    nodes[0].num_edges = edges.size();
    nodes[0].update_bounds(edges, nodes);

    int nodes_used = 1;
    build_nodes_recursive(0, nodes_used);
}

bool BVHTree2::is_inside(const Vector2& p, real dist_min) const {
    //Check multiple directions to minimize chance of hitting an edge

    int num_inside = 0;

    for (int i = 0; i < 4; i++) {
        Vector2 hit_pos;
        Vector2 hit_normal;
        int hit_index;

        nodes[0].ray_cast(p, face_check_dirs[i], edges, nodes, hit_pos, hit_normal, hit_index);
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

    return num_inside >= 2;
}

bool BVHTree2::ray_cast(const Vector2& ray_origin, const Vector2& ray_direction, Vector2& hit_pos, Vector2& hit_normal, int& out_index) const {
    return nodes[0].ray_cast(ray_origin, ray_direction, edges, nodes, hit_pos, hit_normal, out_index);
}

