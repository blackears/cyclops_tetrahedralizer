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


#include "cyclops_tetrahedralizer.h"

#include <random>

using namespace std;
using namespace CyclopsTetra3D;

void Tetrahedron::create_from_points(int v0_idx, int v1_idx, int v2_idx, int v3_idx, const std::vector<Vector3>& points) {
    this->vert_indices[0] = v0_idx;
    this->vert_indices[1] = v1_idx;
    this->vert_indices[2] = v2_idx;
    this->vert_indices[3] = v3_idx;

    this->neighbors[0] = -1;
    this->neighbors[1] = -1;
    this->neighbors[2] = -1;
    this->neighbors[3] = -1;

    this->circumcenter = calc_circumcenter(points);
    this->center = (points[v0_idx] + points[v1_idx] + points[v2_idx] + points[v3_idx]) / 4.0;

    Vector3 p0 = points[v0_idx];
    Vector3 p1 = points[v1_idx];
    Vector3 p2 = points[v2_idx];
    Vector3 p3 = points[v3_idx];

    Plane test_plane(p0, p1, p2);
    if (test_plane.distance_to_plane(p3) > 0.0) {
        //Swap two vertices to change winding
        this->vert_indices[0] = v1_idx;
        this->vert_indices[1] = v0_idx;
    }

    //Should all be facing outside
    for (int i = 0; i < 4; i++) {
        face_planes[i] = Plane(points[face_vert_indices[i][0]], points[face_vert_indices[i][1]], points[face_vert_indices[i][2]]);
        // boundary_face[i] = vert_indices[face_vert_indices[i][0]] < num_boundary_points &&
        //                    vert_indices[face_vert_indices[i][1]] < num_boundary_points &&
        //                    vert_indices[face_vert_indices[i][2]] < num_boundary_points;
    }

    valid = true;
}

Vector3 Tetrahedron::calc_circumcenter(const std::vector<Vector3>& points) const {
    //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter

    //From Matthias Muller
    //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68
    Vector3 p0 = points[vert_indices[0]];
    Vector3 p1 = points[vert_indices[1]];
    Vector3 p2 = points[vert_indices[2]];
    Vector3 p3 = points[vert_indices[3]];

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

bool Tetrahedron::contains_point(const Vector3& p, const std::vector<Vector3>& points) const {
    for (int i = 0; i < 4; i++) {
        real dist = face_planes[i].distance_to_plane(p);
        if (dist <= 0.0) {
            return false;
        }
    }
    return true;
}

int Tetrahedron::find_adjacent_tetrahedron(const Vector3& dir, const std::vector<Vector3>& points) const {
    //constexpr int[][] tet_faces = {{2,1,0}, {0,1,3}, {1,2,3}, {2,0,3}};

    real best_dist = std::numeric_limits<real>::infinity();
    int best_face = -1;
    for (int i = 0; i < 4; i++) {
        Vector3 p_intersect;
        if (face_planes[i].intersect_ray(center, dir, p_intersect)) {
            Vector3 offset = p_intersect - center;
            real dist = offset.dot(dir);
            if (dist > 0.0 && dist < best_dist) {
                best_dist = dist;
                best_face = i;
            }
        }
    }

    return neighbors[best_face];
}

//return value on [0 - 1] where 1 is a perfect tetrahedron
real Tetrahedron::quality(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) const {
    Vector3 d0 = p1 - p0;
    Vector3 d1 = p2 - p0;
    Vector3 d2 = p3 - p0;
    Vector3 d3 = p2 - p1;
    Vector3 d4 = p3 - p2;
    Vector3 d5 = p1 - p3;

    real s0 = d0.magnitude();
    real s1 = d1.magnitude();
    real s2 = d2.magnitude();
    real s3 = d3.magnitude();
    real s4 = d4.magnitude();
    real s5 = d5.magnitude();

    real ms = (s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4 + s5 * s5) / 6.0;
    real rms = sqrt(ms);

    real s = 12.0 / sqrt(2.0);

    real vol = d0.dot(d1.cross(d2)) / 6.0;
    return s * vol / (rms * rms * rms);
}

const Vector3 BVHTree::face_check_dirs[] = {
    Vector3(1, 0, 0),
    Vector3(-1, 0, 0),
    Vector3(0, 1, 0),
    Vector3(0, -1, 0),
    Vector3(0, 0, 1),
    Vector3(0, 0, -1)
};


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
            swap(triangles[i], triangles[j]);
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


void CyclopsTetrahedralizer::create_tetrahedrons(const std::vector<Vector3>& points, 
    const std::vector<int>& indices, 
    float resolution,
    float quality_threshold) {

    //Create BVH from input triangles
    BVHTree bvh_tree;
    bvh_tree.build_from_triangles(points, indices);

    std::vector<Vector3> tess_points;
    tess_points.reserve(points.size() + 4);

    //Add jitter to points to avoid degenerate cases
    std::random_device r;
    std::default_random_engine rng_eng(r());
    std::uniform_real_distribution<real> rand_eps(-1e-5, 1e-5);

    for (const Vector3& p : points) {
        Vector3 jit_p = p + Vector3(rand_eps(rng_eng), rand_eps(rng_eng), rand_eps(rng_eng));
        tess_points.push_back(jit_p);
    }

    //Find bounding box
    Vector3 bb_min = tess_points[0];
    Vector3 bb_max = tess_points[0];
    for (int i = 1; i < tess_points.size(); i++) {
        bb_min = bb_min.min(tess_points[i]);
        bb_max = bb_max.max(tess_points[i]);
    }

    //Add extra points for interior of mesh
    Vector3 bb_size = bb_max - bb_min;
    float max_dim = std::max(bb_size.x, std::max(bb_size.y, bb_size.z));

    if (resolution > 0) {
        int h = max_dim / resolution;

        for (int xi = 0; xi <= int(bb_size.x / h); xi++) {
            float x = bb_min.x + xi * h + rand_eps(rng_eng);
            for (int yi = 0; yi <= int(bb_size.y / h); yi++) {
                float y = bb_min.y + yi * h + rand_eps(rng_eng);
                for (int zi = 0; zi <= int(bb_size.z / h); zi++) {
                    float z = bb_min.z + zi * h + rand_eps(rng_eng);
                    Vector3 p = Vector3(x, y, z);

                    if (bvh_tree.is_inside(p)) {
                        tess_points.push_back(p);
                    }
                }
            }
        }
    }

    //Find bounding tetrahedron
    Vector3 bb_center = (bb_min + bb_max) / 2.0;
    Vector3 btet_v0 = bb_min;
    Vector3 btet_v1 = bb_min + Vector3(bb_size.x * 3.0, 0, 0);
    Vector3 btet_v2 = bb_min + Vector3(0, bb_size.y * 3.0, 0);
    Vector3 btet_v3 = bb_min + Vector3(0, 0, bb_size.z * 3.0);
    //Add margin
    btet_v0 += (btet_v0 - bb_center) * 0.1;
    btet_v1 += (btet_v1 - bb_center) * 0.1;
    btet_v2 += (btet_v2 - bb_center) * 0.1;
    btet_v3 += (btet_v3 - bb_center) * 0.1;

    tess_points.push_back(btet_v0);
    tess_points.push_back(btet_v1);
    tess_points.push_back(btet_v2);
    tess_points.push_back(btet_v3);

    //Create tetrahedrons
    //create_tetrahedrons_internal(tess_points, bvh_tree, quality_threshold);

    std::vector<Tetrahedron> tetrahedrons;
    Tetrahedron tetra;
    tetra.create_from_points(
        int(tess_points.size() - 4),
        int(tess_points.size() - 3),
        int(tess_points.size() - 2),
        int(tess_points.size() - 1),
        points);
    tetrahedrons.push_back(tetra);

    create_tetrahedrons_iter(tetrahedrons, tess_points);

    //Remove exterior tetrahedrons
    for (int i = 0; i < tetrahedrons.size(); i++) {
        Tetrahedron& tet = tetrahedrons[i];
        if (tet.valid) {
            if (!bvh_tree.is_inside(tet.center, 1e-3))
            {
                tet.valid = false;
            }
        }
    }

}

void CyclopsTetrahedralizer::create_tetrahedrons_iter(std::vector<Tetrahedron>& tetrahedrons, const std::vector<Vector3>& points) {
    for (int i = 0; i < points.size() - 4; i++) {
        Vector3 p = points[i];

        int tet_idx = 0;

        //Walk toward containing tetrahedron
        while (tet_idx != -1) {
            Tetrahedron& tet = tetrahedrons[tet_idx];
            if (tet.contains_point(p, points)) {
                break;
            }

            tet_idx = tet.find_adjacent_tetrahedron(p - tet.center, points);
        }

        if (tet_idx == -1) {
            //Could not find containing tetrahedron
            continue;
        }

        //Find tetrahedra with circumspheres containing point
//        bad_tet_candidates.push_back(tet_idx);

        Tetrahedron& tet = tetrahedrons[tet_idx];
        tet.valid = false;

        std::vector<int> bad_tet_indices;

        std::vector<std::tuple<int, int>> bad_tet_candidates;
        for (int i = 0; i < 4; ++i) {
            bad_tet_candidates.push_back(std::make_tuple(tet_idx, i));
        }

        std::vector<std::tuple<int, int>> outer_faces;

        while (!bad_tet_candidates.empty())
        {
            auto [current_tet_idx, face_idx] = bad_tet_candidates.back();
            bad_tet_candidates.pop_back();

            Tetrahedron& current_tet = tetrahedrons[current_tet_idx];

            int neighbor_idx = current_tet.neighbors[face_idx];
            if (neighbor_idx == -1) {
                outer_faces.push_back(std::make_tuple(current_tet_idx, face_idx));
            }
            else {
                Tetrahedron& neighbor_tet = tetrahedrons[neighbor_idx];
                if (!neighbor_tet.valid) {
                    continue;
                }

                if (neighbor_tet.point_in_circumsphere(p, points)) {
                    bad_tet_indices.push_back(neighbor_idx);
                    neighbor_tet.valid = false;

                    for (int i = 0; i < 4; ++i) {
                        bad_tet_candidates.push_back(std::make_tuple(neighbor_idx, i));
                    }
                }
                else {
                    outer_faces.push_back(std::make_tuple(current_tet_idx, face_idx));
                }
            }
        }

        //Rebuild cavity with new tetrahedrons
        std::vector<int> new_tet_indices;
        for (auto [bad_tet_idx, face_idx] : outer_faces) {
            Tetrahedron& bad_tet = tetrahedrons[bad_tet_idx];

            int neighbor_tet_idx = bad_tet.neighbors[face_idx];
            //int v0_idx = 

            int vert_indices[4];
            vert_indices[0] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][0]];
            vert_indices[1] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][1]];
            vert_indices[2] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][2]];
            vert_indices[3] = i;

            Tetrahedron new_tet;
            new_tet.create_from_points(
                vert_indices[0],
                vert_indices[1],
                vert_indices[2],
                vert_indices[3],
                points);
            tetrahedrons.push_back(new_tet);
            int new_tet_idx = tetrahedrons.size() - 1;

            //Update neighbor links to exterior tetrahedrons
            new_tet.neighbors[0] = neighbor_tet_idx;
            //Find face with same vertices with reverse winding
            if (neighbor_tet_idx != -1) {
                Tetrahedron& neighbor_tet = tetrahedrons[neighbor_tet_idx];
                neighbor_tet.neighbors[neighbor_tet.find_face(vert_indices[0], vert_indices[2], vert_indices[1])] = new_tet_idx;
            }

            //Check other cavity filling tetrahedrons for shared faces
            for (int j = 0; j < new_tet_indices.size(); j++) {
                int other_tet_idx = new_tet_indices[j];
                Tetrahedron& other_tet = tetrahedrons[other_tet_idx];

                //Check for shared face
                int shared_count = 0;
                for (int face_idx = 1; face_idx < 4; face_idx++) {
                    int other_face_idx = other_tet.find_face(new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][0]], 
                        new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][1]], 
                        new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][2]]);

                    if (other_face_idx != -1) {
                        //Shared face
                        new_tet.neighbors[face_idx] = other_tet_idx;
                        other_tet.neighbors[other_face_idx] = new_tet_idx;
                    }
                }
            }

            new_tet_indices.push_back(new_tet_idx);
        }

    }
}
