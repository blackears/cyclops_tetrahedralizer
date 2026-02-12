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
#include <algorithm>
#include <fstream>
#include <set>
#include <cassert>
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

    Vector3 p0 = points[v0_idx];
    Vector3 p1 = points[v1_idx];
    Vector3 p2 = points[v2_idx];
    Vector3 p3 = points[v3_idx];

    tet.circumcenter = Math::tetrahedron_circumcenter(p0, p1, p2, p3);
    tet.circumsphere_radius_squared = (tet.circumcenter - p0).magnitude_squared();
    tet.center = (p0 + p1 + p2 + p3) / 4.0;

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

int Tetrahedron::find_adjacent_tetrahedron(const Vector3& dir) const {
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

    if (best_face == -1)
        return -1;

    return neighbors[best_face];
}

int Tetrahedron::step_toward_point_adjacent_tetrahedron(const Vector3& p, real epsilon) const {
    Vector3 p_offset = p - center;

    real best_dist_sq = std::numeric_limits<real>::infinity();
    int best_face = -1;
    for (int i = 0; i < 4; i++) {
        Vector3 f_intersect;
        if (face_planes[i].intersect_ray(center, p_offset, f_intersect)) {
            Vector3 isect_offset = (f_intersect - center);
            if (isect_offset.magnitude_squared() < epsilon * epsilon)
                return -1;

            if (isect_offset.dot(p_offset) <= 0.0)
                continue;

            real dist_sq = isect_offset.magnitude_squared();
            if (dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best_face = i;
            }
        }
    }

    if (best_face == -1 || best_dist_sq > p_offset.magnitude_squared())
        return -1;

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
    float subdivisions) {

    //Create BVH from input triangles
    BVHTree3 bvh_tree;
    bvh_tree.build_from_triangles(points, indices);

    tess_points.clear();
    tess_points.reserve(points.size() + 4);

    //Add jitter to points to avoid degenerate cases
    std::random_device rd;
    std::default_random_engine rng_eng(rd());
    rng_eng.seed(0);
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

    if (subdivisions > 0) {
        real max_dim = std::max(bb_size.x, std::max(bb_size.y, bb_size.z));
        real cube_side_len = max_dim / subdivisions;

        int steps_x = ceil(bb_size.x / cube_side_len);
        int steps_y = ceil(bb_size.y / cube_side_len);
        int steps_z = ceil(bb_size.z / cube_side_len);

        Vector3 grid_size(steps_x * cube_side_len, steps_y * cube_side_len, steps_z * cube_side_len);

        for (int k = 0; k < steps_z; ++k) {
            for (int j = 0; j < steps_y; ++j) {
                for (int i = 0; i < steps_x; ++i) {
                    Vector3 p = Vector3(i, j, k) * cube_side_len - grid_size / 2 + bb_min + bb_size / 2;
                    p += Vector3(rand_eps(rng_eng), rand_eps(rng_eng), rand_eps(rng_eng));
                    
                    if (bvh_tree.is_inside(p, false)) {
                        tess_points.push_back(p);
                    }

                }
            }
        }
    }

    std::mt19937 m_eng(rd());
    std::shuffle(tess_points.begin(), tess_points.end(), m_eng);

    //Find bounding tetrahedron
    Vector3 bb_center = (bb_min + bb_max) / 2.0;
    Vector3 btet_v0 = bb_min;
    Vector3 btet_v1 = bb_min + Vector3(bb_size.x * 3.0, 0, 0);
    Vector3 btet_v2 = bb_min + Vector3(0, bb_size.y * 3.0, 0);
    Vector3 btet_v3 = bb_min + Vector3(0, 0, bb_size.z * 3.0);
    //Add margin
    btet_v0 += (btet_v0 - bb_center) * 0.2;
    btet_v1 += (btet_v1 - bb_center) * 0.2;
    btet_v2 += (btet_v2 - bb_center) * 0.2;
    btet_v3 += (btet_v3 - bb_center) * 0.2;

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

    create_tetrahedrons_iter(tess_points);

    ////////////////
//    bool inside = bvh_tree.is_inside(Vector3(0, .216, 0), 1e-3);
    //bool inside = bvh_tree.is_inside(tetrahedra[13394].center, 1e-3);
    //int j = 9;
    //bool inside2 = bvh_tree.is_inside(tetrahedra[13394].center, 1e-3);
    ////////////////

    //Remove exterior tetrahedrons
    for (int i = 0; i < tetrahedra.size(); i++) {
        Tetrahedron& tet = tetrahedra[i];
        if (tet.valid) {
            if (!bvh_tree.is_inside(tet.center, true))
            {
                tet.valid = false;
            }
        }
    }

}

void CyclopsTetrahedralizer::create_tetrahedrons_iter(const std::vector<Vector3>& points) {
    //Last 4 points are bounding tetrahedron
    for (int p_idx = 0; p_idx < points.size() - 4; p_idx++) {
        Vector3 p = points[p_idx];

        int tet_idx = 0;

        //Skip forward to first valid tetrahedron
        while (tet_idx != -1) {
            Tetrahedron& tri = tetrahedra[tet_idx];
            if (tri.valid)
                break;

            tet_idx++;
        }

        //Walk toward containing tetrahedron
        while (tet_idx != -1) {
            Tetrahedron& tet = tetrahedra[tet_idx];
            int next_tet_idx = tet.step_toward_point_adjacent_tetrahedron(p);
            if (next_tet_idx == -1)
                break;

            tet_idx = next_tet_idx;
        }

        if (tet_idx == -1) {
            //Could not find containing tetrahedron
            continue;
        }

        //Find tets which have circumcenters that include point
        std::vector<int> tets_to_scan;
        std::vector<int> tets_to_replace;
        std::set<int> tets_viewed;
        tets_to_scan.push_back(tet_idx);
        tets_to_replace.push_back(tet_idx);
        tets_viewed.emplace(tet_idx);

        while (!tets_to_scan.empty()) {
            int cur_tet_idx = tets_to_scan.back();
            Tetrahedron& current_tet = tetrahedra[cur_tet_idx];
            tets_to_scan.pop_back();

            for (int i = 0; i < 4; ++i) {
                int neighbor_tet_idx = current_tet.neighbors[i];
                
                if (neighbor_tet_idx == -1 || tets_viewed.find(neighbor_tet_idx) != tets_viewed.end())
                    continue;

                tets_viewed.emplace(neighbor_tet_idx);
                
                Tetrahedron& neighbor_tet = tetrahedra[neighbor_tet_idx];
                if ((neighbor_tet.circumcenter - p).magnitude_squared() < neighbor_tet.circumsphere_radius_squared) {
                    tets_to_replace.push_back(neighbor_tet_idx);
                    tets_to_scan.push_back(neighbor_tet_idx);
                }
            }
        }

        std::vector<std::tuple<int, int>> boundary_faces;
        while (true) {
            //Find cavity bounds
            boundary_faces.clear();
            for (int cur_tet_idx : tets_to_replace) {
                Tetrahedron& current_tet = tetrahedra[cur_tet_idx];

                for (int i = 0; i < 4; ++i) {
                    int neighbor_tet_idx = current_tet.neighbors[i];

                    if (neighbor_tet_idx == -1 || std::find(tets_to_replace.begin(), tets_to_replace.end(), neighbor_tet_idx) == tets_to_replace.end()) {
                        boundary_faces.push_back(std::tuple<int, int>(cur_tet_idx, i));
                    }
                }
            }

            //Boundary faces should be a convex shape, but due to round off errors may have concavity
            //Remove tets that create a concave boundary
            std::set<int> tets_violating;
            for (auto [cur_tet_idx, boundary_face_idx] : boundary_faces) {
                Tetrahedron& cur_tet = tetrahedra[cur_tet_idx];

                int neighbor_tet_idx = cur_tet.neighbors[boundary_face_idx];

                int vi_0 = cur_tet.vert_indices[Tetrahedron::face_vert_indices[boundary_face_idx][0]];
                int vi_1 = cur_tet.vert_indices[Tetrahedron::face_vert_indices[boundary_face_idx][1]];
                int vi_2 = cur_tet.vert_indices[Tetrahedron::face_vert_indices[boundary_face_idx][2]];
                int vi_3 = p_idx;

                const Vector3& p0 = tess_points[vi_0];
                const Vector3& p1 = tess_points[vi_1];
                const Vector3& p2 = tess_points[vi_2];
                const Vector3& p3 = tess_points[vi_3];

                real volume_x2 = Math::det(p0 - p3, p1 - p3, p2 - p3);
                if (volume_x2 <= 0) {
                    tets_violating.emplace(cur_tet_idx);
                }
            }

            if (tets_violating.size() == 0)
                break;

            //save_file_obj("concavity.obj");

            //Negative volumes indicate a concavity
            //Remove all tets that had negative volume
            tets_to_replace.erase(std::remove_if(tets_to_replace.begin(), tets_to_replace.end(),
                    [&](int x) { return tets_violating.count(x) > 0; }
                ), tets_to_replace.end());
        }

        if (tets_to_replace.empty())
            continue;

        //Mark invalid
        for (int tet_idx : tets_to_replace) {
            Tetrahedron& current_tet = tetrahedra[tet_idx];
            current_tet.valid = false;
        }

        //Add new tets
        std::vector<int> tets_added;
        for (auto [bad_tet_idx, bad_face_idx] : boundary_faces) {
            Tetrahedron& bad_tet = tetrahedra[bad_tet_idx];

            int neighbor_tet_idx = bad_tet.neighbors[bad_face_idx];

            int vi_0 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][0]];
            int vi_1 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][1]];
            int vi_2 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][2]];
            int vi_3 = p_idx;

            int new_tet_idx = tetrahedra.size();
            tets_added.push_back(new_tet_idx);

            //First face should match winding of outer face of cur_tet
            tetrahedra.push_back(Tetrahedron::create_from_points(vi_0, vi_1, vi_2, vi_3, points));
            Tetrahedron& new_tet = tetrahedra[new_tet_idx];

            //Should have already filtered out tets with negative volume
            assert(new_tet.volume_times_2(points) > 0);

            if (neighbor_tet_idx != -1) {
                Tetrahedron& neighbor_tet = tetrahedra[neighbor_tet_idx];

                int neighbor_face_idx = neighbor_tet.find_face(vi_0, vi_2, vi_1);
                neighbor_tet.neighbors[neighbor_face_idx] = new_tet_idx;
                new_tet.neighbors[0] = neighbor_tet_idx;
            }
        }

        //Set neighbors of added tets
        for (int i = 0; i < tets_added.size() - 1; ++i) {
            int tet_0_idx = tets_added[i];
            Tetrahedron& tet_0 = tetrahedra[tet_0_idx];

            for (int j = i + 1; j < tets_added.size(); ++j) {
                int tet_1_idx = tets_added[j];
                Tetrahedron& tet_1 = tetrahedra[tet_1_idx];

                //Face 0 on all added tets faces boundary
                for (int k = 1; k < 4; ++k) {
                    int vi_0 = tet_0.vert_indices[Tetrahedron::face_vert_indices[k][0]];
                    int vi_1 = tet_0.vert_indices[Tetrahedron::face_vert_indices[k][1]];
                    int vi_2 = tet_0.vert_indices[Tetrahedron::face_vert_indices[k][2]];

                    int neighbor_face = tet_1.find_face(vi_0, vi_2, vi_1);
                    if (neighbor_face != -1) {
                        tet_0.neighbors[k] = tet_1_idx;
                        tet_1.neighbors[neighbor_face] = tet_0_idx;
                        //Tets will only match on one face
                        break;
                    }
                }
            }
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

void CyclopsTetrahedralizer::save_file_line_segments_obj(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# Cyclops Tetrahedralizer" << std::endl;
    file << "# https://github.com/blackears/cyclops_tetrahedralizer" << std::endl;
    int p_idx = 0;
    for (const auto& p : tess_points) {
        file << "v " << p.x << " " << p.y << " " << p.z << " \t#" << p_idx++ + 1 << std::endl;
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
                    file << "l " << vi0 + 1 << " " << vi1 + 1 << std::endl;
                }
            }
        }
    }

    file.close();
}

void CyclopsTetrahedralizer::save_file_obj(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# Cyclops Tetrahedralizer" << std::endl;
    file << "# https://github.com/blackears/cyclops_tetrahedralizer" << std::endl;

    int p_idx = 0;
    for (const auto& p : tess_points) {
        file << "v " << p.x << " " << p.y << " " << p.z << " \t#" << p_idx++ + 1 << std::endl;
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

    for (int tet_idx = 0; tet_idx < tetrahedra.size(); ++tet_idx) {
        const Tetrahedron& tet = tetrahedra[tet_idx];
        if (!tet.valid)
            continue;

        file << "o tet_" << tet_idx << std::endl;

        const std::array<int, 4> vi = tet.get_vert_indices();

        file << "f " << vi[0] + 1 << "/1 " << vi[1] + 1 << "/2 " << vi[2] + 1 << "/3" << std::endl;
        file << "f " << vi[1] + 1 << "/4 " << vi[0] + 1 << "/5 " << vi[3] + 1 << "/6" << std::endl;
        file << "f " << vi[2] + 1 << "/7 " << vi[3] + 1 << "/8 " << vi[0] + 1 << "/9" << std::endl;
        file << "f " << vi[3] + 1 << "/10 " << vi[2] + 1 << "/11 " << vi[1] + 1 << "/12" << std::endl;
    }

    file.close();
}

void CyclopsTetrahedralizer::dump_outer_faces(const std::vector<std::tuple<int, int>>& outer_faces) {
    for (auto [bad_tet_idx, bad_face_idx] : outer_faces) {
        Tetrahedron& bad_tet = tetrahedra[bad_tet_idx];

        int vi_0 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][0]];
        int vi_1 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][1]];
        int vi_2 = bad_tet.vert_indices[Tetrahedron::face_vert_indices[bad_face_idx][2]];

        std::cout << tess_points[vi_0] << ", " << tess_points[vi_1] << ", " << tess_points[vi_2] << ", #tet_idx " << bad_tet_idx << "  face_idx " << bad_face_idx << std::endl;
    }

}
