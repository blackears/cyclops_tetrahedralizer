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


#include "cyclops_tessellate.h"

#include <random>

using namespace std;
using namespace CyclopsTessellate3D;

//https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py

template<typename T>
Vector3<T> Tetrahedron<T>::calc_circum_center(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const {
    //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter

    //From Matthias Muller
    //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68

    Vector3<T> b = p1 - p0;
    Vector3<T> c = p2 - p0;
    Vector3<T> d = p3 - p0;
 
    T det = 2.0 * (b.x * (c.y * d.z - c.z * d.y) 
        - b.y * (c.x * d.z - c.z * d.x) 
        + b.z * (c.x * d.y - c.y * d.x));

    if (det == 0.0) {
        return p0;
    }
    else {
        Vector3<T> v = c.cross(d) * b.dot(b) + d.cross(b) * c.dot(c) + b.cross(c) * d.dot(d);
        v /= det;
        return p0 + v;
    }
}

//return value on [0 - 1] where 1 is a perfect tetrahedron
template<typename T>
T Tetrahedron<T>::quality(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const {
    Vector3<T> d0 = p1 - p0;
    Vector3<T> d1 = p2 - p0;
    Vector3<T> d2 = p3 - p0;
    Vector3<T> d3 = p2 - p1;
    Vector3<T> d4 = p3 - p2;
    Vector3<T> d5 = p1 - p3;

    T s0 = d0.magnitude();
    T s1 = d1.magnitude();
    T s2 = d2.magnitude();
    T s3 = d3.magnitude();
    T s4 = d4.magnitude();
    T s5 = d5.magnitude();

    T ms = (s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4 + s5 * s5) / 6.0;
    T rms = sqrt(ms);

    T s = 12.0 / sqrt(2.0);

    T vol = d0.dot(d1.cross(d2)) / 6.0;
    return s * vol / (rms * rms * rms);
}


template<typename T>
void BVHTree<T>::build_nodes_recursive(int root_node_idx, int& nodes_used) {
    BVHTreeNode<T>& root_node = nodes[root_node_idx];

    Vector3<T> bb_size = root_node.bounds.size();
    int axis;
    T midpoint
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
    int i = root_node.first_primitive_offset;
    int j = i + root_node.num_triangles - 1;
    while (i <= j) {
        BVHTreeTriangle<T>& tri = triangles[i];
        if (tri.centroid[axis] < midpoint) {
            i++;
        } else {
            swap(triangles[i], triangles[j]);
            j--;
        }
    }
    
    int left_count = i - root_node.first_triangle_offset;
    if (left_count == 0 || left_count == root_node.num_triangles) {
        return; //Cannot split further
    }

    //Next two nodes in node list with be children of this node
    nodes[nodes_used].first_primitive_offset = root_node.first_triangle_offset;
    nodes[nodes_used].num_primitives = left_count;
    nodes[nodes_used].update_bounds(triangles, nodes);

    nodes[nodes_used + 1].first_primitive_offset = i;
    nodes[nodes_used + 1].num_primitives = root_node.num_triangles - left_count;
    nodes[nodes_used + 1].update_bounds(triangles, nodes);
    
    root_node.left_child_idx = nodes_used;
    root_node.num_triangles = 0;
    nodes_used += 2;

    //Recurse
    build_nodes_recursive(root_node.child_left_idx, nodes_used);
    build_nodes_recursive(root_node.child_left_idx + 1, nodes_used);
}

template<typename T>
void BVHTree<T>::build_from_triangles(const std::vector<Vector3<T>>& points, const std::vector<int>& indices) {
    nodes.clear();
    triangles.clear();

    nodes.resize(indices.size() / 3 * 2 - 1);
    triangles.reserve(indices.size() / 3);

    for (int i = 0; i < indices.size() / 3; i++) {
        int idx0 = indices[i * 3 + 0];
        int idx1 = indices[i * 3 + 1];
        int idx2 = indices[i * 3 + 2];

        Vector3<T> v0 = points[idx0];
        Vector3<T> v1 = points[idx1];
        Vector3<T> v2 = points[idx2];

        BVHTreeTriangle<T> tri = BVHTreeTriangle<T>(v0, v1, v2);
        triangles.push_back(tri);
    }

    nodes[0].left_child_idx = 0;
    nodes[0].right_child_idx = 0;
    nodes[0].first_primitive_offset = 0;
    nodes[0].num_primitives = triangles.size();
    nodes[0].update_bounds();

    int nodes_used = 1;
    build_nodes_recursive(0, nodes_used);
}

template<typename T>
bool BVHTree<T>::is_inside(Vector3<T> p, T dist_min) const {
    //Check multiple directions to minimize chance of hitting an edge
    const Vector3[T] dirs = {
        Vector3<T>(1,0,0),
        Vector3<T>(-1,0,0),
        Vector3<T>(0,1,0),
        Vector3<T>(0,-1,0),
        Vector3<T>(0,0,1),
        Vector3<T>(0,0,-1)
    };

    int num_inside = 0;

    for (int i = 0; i < 6; i++) {
        Vector3<T> hit_pos;
        Vector3<T> hit_normal;
        int hit_index;

        nodes[0].ray_cast(p, dirs[i], 1e10, hit_pos, hit_normal, hit_index);
        T hit_distance = (hit_pos - p).magnitude();
        
        //Check hit is valid
        if (hit_index >= 0) {
            if (hit_normal.dot(dirs[i]) > 0) {
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

template<typename T>
bool BVHTree<T>::ray_cast(Vector3<T> ray_origin, Vector3<T> ray_direction, Vector3<T> &hit_pos, Vector3<T> &hit_normal, int &out_index) const {
    return nodes[0].ray_cast(ray_origin, ray_direction, hit_pos, hit_normal, out_index);
}


template<typename T>
void create_tetrahedron_ids(const std::vector<Vector3<T>>& points, BVHTree<T>& bvh_tree, float quality_threshold) {

}

template<typename T>
void CyclopsTetrahedralizer<T>::create_tetrahedrons(const std::vector<Vector3<T>>& points, 
    const std::vector<int>& indices, 
    float resolution,
    float quality_threshold) {

    //Create BVH from input triangles
    BVHTree<T> bvh_tree;
    bvh_tree.build_from_triangles(points, indices);


    std::vector<Vector3<T>> tess_points;

    //Add jitter to points to avoid degenerate cases
    std::random_device r;
    std::default_random_engine rng_eng(r());
    std::uniform_real_distribution<T> rand_eps(-1e-5, 1e-5);

    tess_points.reserve(points.size());
    for (const Vector3<T>& p : points) {
        Vector3<T> jit_p = p + Vector3<T>(rand_eps(rng_eng), rand_eps(rng_eng), rand_eps(rng_eng));
        tess_points.push_back(jit_p);
    }

    //Find bounding box
    Vector3<T> bb_min = tess_points[0];
    Vector3<T> bb_max = tess_points[0];
    for (const Vector3<T>& p : tess_points) {
        bb_min = bb_min.min(p);
        bb_max = bb_max.max(p);
    }

    //Add extra points for interior of mesh
    Vector3<T> bb_size = bb_max - bb_min;
    float max_dim = std::max(bb_size.x, std::max(bb_size.y, bb_size.z));

    if (resolution > 0) {
        int h = max_dim / resolution;

        for (int xi = 0; xi <= int(bb_size.x / h); xi++) {
            float x = bb_min.x + xi * h + rand_eps(rng_eng);
            for (int yi = 0; yi <= int(bb_size.y / h); yi++) {
                float y = bb_min.y + yi * h + rand_eps(rng_eng);
                for (int zi = 0; zi <= int(bb_size.z / h); zi++) {
                    float z = bb_min.z + zi * h + rand_eps(rng_eng);
                    Vector3<T> p = Vector3<T>(x, y, z);

                    if (bvh_tree.is_inside(p)) {
                        tess_points.push_back(p);
                    }
                }
            }
        }
    }

    //Find bounding tetrahedron
    Vector3<T> bb_center = (bb_min + bb_max) / 2.0;
    //Vector3<T> bb_size = bb_max - bb_min;
    Vector3<T> btet_v0 = bb_min;
    Vector3<T> btet_v1 = bb_min + Vector3<T>(bb_size.x * 3.0, 0, 0);
    Vector3<T> btet_v2 = bb_min + Vector3<T>(0, bb_size.y * 3.0, 0);
    Vector3<T> btet_v3 = bb_min + Vector3<T>(0, 0, bb_size.z * 3.0);
    //Add margin
    btet_v0 += (btet_v0 - bb_center) * 0.1;
    btet_v1 += (btet_v1 - bb_center) * 0.1;
    btet_v2 += (btet_v2 - bb_center) * 0.1;
    btet_v3 += (btet_v3 - bb_center) * 0.1;

    //Create tetrahedrons
    create_tetrahedron_ids(tess_points, bvh_tree, quality_threshold);

}