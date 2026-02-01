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
#include <fstream>
#include <set>
#include "bvh_tree3.h"

using namespace CyclopsTetra3D;

Tetrahedron Tetrahedron::create_from_points(int v0_idx, int v1_idx, int v2_idx, int v3_idx, const std::vector<Vector3>& points) {
    Tetrahedron tet;

    tet.vert_indices[0] = v0_idx;
    tet.vert_indices[1] = v1_idx;
    tet.vert_indices[2] = v2_idx;
    tet.vert_indices[3] = v3_idx;

    tet.neighbors[0] = -1;
    tet.neighbors[1] = -1;
    tet.neighbors[2] = -1;
    tet.neighbors[3] = -1;

    tet.circumcenter = Math::tetrahedron_circumcenter(points[v0_idx], points[v1_idx], points[v2_idx], points[v3_idx]);
    tet.circumcircle_radius_squared = (tet.circumcenter - points[v0_idx]).magnitude_squared();
    tet.center = (points[v0_idx] + points[v1_idx] + points[v2_idx] + points[v3_idx]) / 4.0;

    Vector3 p0 = points[v0_idx];
    Vector3 p1 = points[v1_idx];
    Vector3 p2 = points[v2_idx];
    Vector3 p3 = points[v3_idx];

    Plane test_plane = Plane::create(p0, p1, p2);
    if (test_plane.distance_to_plane(p3) > 0.0) {
        //Swap two vertices to change winding
        tet.vert_indices[0] = v1_idx;
        tet.vert_indices[1] = v0_idx;
    }

    //Should all be facing outside
    for (int i = 0; i < 4; i++) {
        const Vector3& pl_p0 = points[tet.vert_indices[face_vert_indices[i][0]]];
        const Vector3& pl_p1 = points[tet.vert_indices[face_vert_indices[i][1]]];
        const Vector3& pl_p2 = points[tet.vert_indices[face_vert_indices[i][2]]];
        tet.face_planes[i] = Plane::create(pl_p0, pl_p1, pl_p2);
    }

    tet.valid = true;
    return tet;
}

bool Tetrahedron::contains_point(const Vector3& p, const std::vector<Vector3>& points) const {
    return Math::tetrahedron_contains_point(p, points[vert_indices[0]], points[vert_indices[1]], points[vert_indices[2]], points[vert_indices[3]]);
}

int Tetrahedron::find_adjacent_tetrahedron(const Vector3& dir, const std::vector<Vector3>& points) const {
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


void CyclopsTetrahedralizer::create_tetrahedrons(const std::vector<Vector3>& points, 
    const std::vector<int>& indices, 
    float resolution) {

    //Create BVH from input triangles
    BVHTree3 bvh_tree;
    bvh_tree.build_from_triangles(points, indices);

    tess_points.clear();
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

    if (resolution > 0) {
        float max_dim = std::max(bb_size.x, std::max(bb_size.y, bb_size.z));
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

    //Create bounding tetrahedron - reverse winding
    tetrahedra.push_back(Tetrahedron::create_from_points(
        int(tess_points.size() - 4),
        int(tess_points.size() - 2),
        int(tess_points.size() - 3),
        int(tess_points.size() - 1),
        tess_points));

    create_tetrahedrons_iter(tetrahedra, tess_points);

    //Skip exterior removal for now
    return;

    //Remove exterior tetrahedrons
    for (int i = 0; i < tetrahedra.size(); i++) {
        Tetrahedron& tet = tetrahedra[i];
        if (tet.valid) {
            if (!bvh_tree.is_inside(tet.center, 1e-3))
            {
                tet.valid = false;
            }
        }
    }

}

void CyclopsTetrahedralizer::create_tetrahedrons_iter(std::vector<Tetrahedron>& tetrahedrons, const std::vector<Vector3>& points) {
    //Last 4 points are bounding tetrahedron
    for (int i = 0; i < points.size() - 4; i++) {
        Vector3 p = points[i];

        int tet_idx = 0;

        //Skip forward to first valid tetrahedron
        while (tet_idx != -1) {
            Tetrahedron& tri = tetrahedrons[tet_idx];
            if (tri.valid)
                break;

            tet_idx++;
        }

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

            int vert_indices[4];
            //Reverse winding for adjacent tetrahedron
            vert_indices[0] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][0]];
            vert_indices[1] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][2]];
            vert_indices[2] = bad_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][1]];
            vert_indices[3] = i;

            int new_tet_idx = tetrahedrons.size();
            tetrahedrons.push_back(Tetrahedron::create_from_points(
                vert_indices[0],
                vert_indices[1],
                vert_indices[2],
                vert_indices[3],
                points));

            Tetrahedron& new_tet = tetrahedrons[new_tet_idx];

            //Update neighbor links to exterior tetrahedrons
            new_tet.neighbors[0] = neighbor_tet_idx;
            //Find face with same vertices with reverse winding
            if (neighbor_tet_idx != -1) {
                Tetrahedron& neighbor_tet = tetrahedrons[neighbor_tet_idx];
                int face_idx = neighbor_tet.find_face(
                    new_tet.vert_indices[Tetrahedron::face_vert_indices[0][0]], 
                    new_tet.vert_indices[Tetrahedron::face_vert_indices[0][2]], 
                    new_tet.vert_indices[Tetrahedron::face_vert_indices[0][1]]);

                neighbor_tet.neighbors[face_idx] = new_tet_idx;
            }

            //Check other cavity filling tetrahedrons for shared faces
            for (int j = 0; j < new_tet_indices.size(); j++) {
                int other_tet_idx = new_tet_indices[j];
                Tetrahedron& other_tet = tetrahedrons[other_tet_idx];

                //Check for shared face
                for (int face_idx = 1; face_idx < 4; face_idx++) {
                    int vi_0 = new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][0]];
                    int vi_1 = new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][1]];
                    int vi_2 = new_tet.vert_indices[Tetrahedron::face_vert_indices[face_idx][2]];
                    //Reverse winding
                    int other_face_idx = other_tet.find_face(vi_0, vi_2, vi_1);

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

void CyclopsTetrahedralizer::get_mesh(std::vector<Vector3>& out_points, std::vector<int>& out_indices) {
    out_indices.resize(tetrahedra.size() * 12);

    int count = 0;
    for (Tetrahedron& tet : tetrahedra) {
        for (int i = 0; i < 4; ++i) {
            out_indices[count++] = tet.get_vert_indices()[i];
            out_indices[count++] = tet.get_vert_indices()[i];
            out_indices[count++] = tet.get_vert_indices()[i];
            out_indices[count++] = tet.get_vert_indices()[i];

        }
    }
}

void CyclopsTetrahedralizer::save_obj_file_line_segments(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# Cyclops Tetrahedralizer" << std::endl;
    file << "# https://github.com/blackears/cyclops_tetrahedralizer" << std::endl;
    for (const auto& p : tess_points) {
        file << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }

    std::set<Vector2> used_edges;
    int tet_count = 0;
    for (auto& tet : tetrahedra) {
        if (!tet.valid)
            continue;

        for (int i = 0; i <= 2; ++i) {
            for (int j = i + 1; j <= 3; ++j) {
                int vi0 = tet.vert_indices[i];
                int vi1 = tet.vert_indices[j];

                if (used_edges.find(Vector2(vi0, vi1)) == used_edges.end() && used_edges.find(Vector2(vi1, vi0)) == used_edges.end()) {
                    used_edges.insert(Vector2(vi0, vi1));
                    file << "l " << vi0 << " " << vi1 << std::endl;
                }
            }
        }
    }

    file.close();
}

void CyclopsTetrahedralizer::save_obj_file(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# Cyclops Tetrahedralizer" << std::endl;
    file << "# https://github.com/blackears/cyclops_tetrahedralizer" << std::endl;
    for (const auto& p : tess_points) {
        file << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }

    for (auto& tet : tetrahedra) {
        if (!tet.valid)
            continue;

        for (int i = 0; i < 4; ++i) {
            const Vector3& n = tet.face_planes[i].normal;
            file << "vn " << n.x << " " << n.y << " " << n.z << std::endl;
        }
    }

    file << "vt 0 0" << std::endl;
    file << "vt .5 0" << std::endl;
    file << "vt .25 .5" << std::endl;
    file << "vt .5 0" << std::endl;
    file << "vt 1 0" << std::endl;
    file << "vt .25 .5" << std::endl;
    file << "vt 0 .5" << std::endl;
    file << "vt .5 .5" << std::endl;
    file << "vt .25 1" << std::endl;
    file << "vt .5 .5" << std::endl;
    file << "vt 1 .5" << std::endl;
    file << "vt .25 1" << std::endl;

    int tet_count = 0;
    for (auto& tet : tetrahedra) {
        if (!tet.valid)
            continue;

        for (int j = 0; j < 4; ++j) {
            file << "f";

            for (int i = 0; i < 3; ++i) {
                file << " " << tet.vert_indices[Tetrahedron::face_vert_indices[j][i]]
                    << "/" << (j * 3 + i)
                    << "/" << (tet_count * 4 + j);
            }

            file << std::endl;
        }

        tet_count++;
    }

    file.close();
}
