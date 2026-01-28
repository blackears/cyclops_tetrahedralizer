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

#include "delunay_triangulation.h"

#include <random>
#include "bvh_tree2.h"

using namespace CyclopsTetra3D;

void DelunayTriangulator::create_triangles(const std::vector<Vector2>& points, const std::vector<int>& indices, float resolution) {
    //Create BVH from input triangles
    BVHTree2 bvh_tree;
    bvh_tree.build_from_edges(points, indices);

    std::vector<Vector2> tess_points;
    tess_points.reserve(points.size() + 4);

    //Add jitter to points to avoid degenerate cases
    std::random_device r;
    std::default_random_engine rng_eng(r());
    std::uniform_real_distribution<real> rand_eps(-1e-5, 1e-5);

    for (const Vector2& p : points) {
        Vector2 jit_p = p + Vector2(rand_eps(rng_eng), rand_eps(rng_eng));
        tess_points.push_back(jit_p);
    }

    //Find bounding box
    Vector2 bb_min = tess_points[0];
    Vector2 bb_max = tess_points[0];
    for (int i = 1; i < tess_points.size(); i++) {
        bb_min = bb_min.min(tess_points[i]);
        bb_max = bb_max.max(tess_points[i]);
    }

    Vector2 bb_size = bb_max - bb_min;

    if (resolution > 0) {
        float max_dim = std::max(bb_size.x, bb_size.y);
        int h = max_dim / resolution;

        for (int xi = 0; xi <= int(bb_size.x / h); xi++) {
            float x = bb_min.x + xi * h + rand_eps(rng_eng);
            for (int yi = 0; yi <= int(bb_size.y / h); yi++) {
                float y = bb_min.y + yi * h + rand_eps(rng_eng);
                Vector2 p = Vector2(x, y);

                if (bvh_tree.is_inside(p)) {
                    tess_points.push_back(p);
                }
            }
        }
    }

    //Find bounding triangle
    Vector2 bb_center = (bb_min + bb_max) / 2.0;
    Vector2 btri_v0 = bb_min;
    Vector2 btri_v1 = bb_min + Vector2(bb_size.x * 3.0, 0);
    Vector2 btri_v2 = bb_min + Vector2(0, bb_size.y * 3.0);
    //Add margin
    btri_v0 += (btri_v0 - bb_center) * 0.1;
    btri_v1 += (btri_v1 - bb_center) * 0.1;
    btri_v2 += (btri_v2 - bb_center) * 0.1;

    tess_points.push_back(btri_v0);
    tess_points.push_back(btri_v1);
    tess_points.push_back(btri_v2);

    triangles.push_back(DelunayTriangle::create_from_points(
        int(tess_points.size() - 3),
        int(tess_points.size() - 2),
        int(tess_points.size() - 1),
        tess_points));

    create_triangles_iter(tess_points);

    //Remove exterior triangles
    for (int i = 0; i < triangles.size(); i++) {
        DelunayTriangle& tri = triangles[i];
        if (tri.valid) {
            if (!bvh_tree.is_inside(tri.center, 1e-3))
            {
                tri.valid = false;
            }
        }
    }
}

void DelunayTriangulator::create_triangles_iter(const std::vector<Vector2>& points) {
    for (int i = 0; i < points.size() - 3; i++) {
        Vector2 p = points[i];

        int tri_idx = 0;

        while (tri_idx != -1) {
            DelunayTriangle& tri = triangles[tri_idx];
            if (tri.valid)
                break;
            
            tri_idx++;
        }

        //Walk toward containing triangle
        while (tri_idx != -1) {
            DelunayTriangle& tri = triangles[tri_idx];
            if (tri.contains_point(p, points)) {
                break;
            }

            tri_idx = tri.find_adjacent_triangle(p - tri.center, points);
        }

        if (tri_idx == -1) {
            //Could not find containing triangle
            continue;
        }

        //Find triangles with circumspheres containing point
        DelunayTriangle& tri = triangles[tri_idx];
        tri.valid = false;

        std::vector<int> bad_tri_indices;

        std::vector<std::tuple<int, int>> bad_tri_candidates;
        for (int i = 0; i < 3; ++i) {
            bad_tri_candidates.push_back(std::make_tuple(tri_idx, i));
        }

        std::vector<std::tuple<int, int>> outer_edges;

        while (!bad_tri_candidates.empty())
        {
            auto [current_tri_idx, edge_idx] = bad_tri_candidates.back();
            bad_tri_candidates.pop_back();

            DelunayTriangle& current_tri = triangles[current_tri_idx];

            int neighbor_idx = current_tri.neighbors[edge_idx];
            if (neighbor_idx == -1) {
                outer_edges.push_back(std::make_tuple(current_tri_idx, edge_idx));
            }
            else {
                DelunayTriangle& neighbor_tri = triangles[neighbor_idx];
                if (!neighbor_tri.valid) {
                    continue;
                }

                if (neighbor_tri.point_in_circumcircle(p, points)) {
                    bad_tri_indices.push_back(neighbor_idx);
                    neighbor_tri.valid = false;

                    for (int i = 0; i < 3; ++i) {
                        bad_tri_candidates.push_back(std::make_tuple(neighbor_idx, i));
                    }
                }
                else {
                    outer_edges.push_back(std::make_tuple(current_tri_idx, edge_idx));
                }
            }
        }

        //Rebuild cavity with new triangles
        std::vector<int> new_tri_indices;
        for (auto [bad_tri_idx, edge_idx] : outer_edges) {
            DelunayTriangle& bad_tri = triangles[bad_tri_idx];

            int neighbor_tri_idx = bad_tri.neighbors[edge_idx];

            int vert_indices[3];
            vert_indices[0] = bad_tri.vert_indices[DelunayTriangle::edge_vert_indices[edge_idx][0]];
            vert_indices[1] = bad_tri.vert_indices[DelunayTriangle::edge_vert_indices[edge_idx][1]];
            vert_indices[2] = i;

            int new_tri_idx = triangles.size();
            triangles.push_back(DelunayTriangle::create_from_points(
                vert_indices[0],
                vert_indices[1],
                vert_indices[2],
                points));

            DelunayTriangle& new_tri = triangles[new_tri_idx];

            //Update neighbor links to exterior triangles
            new_tri.neighbors[0] = neighbor_tri_idx;
            //Find face with same vertices with reverse winding
            if (neighbor_tri_idx != -1) {
                DelunayTriangle& neighbor_tri = triangles[neighbor_tri_idx];
                neighbor_tri.neighbors[neighbor_tri.find_edge(vert_indices[1], vert_indices[0])] = new_tri_idx;
            }

            //Check other cavity filling triangles for shared faces
            for (int j = 0; j < new_tri_indices.size(); j++) {
                int other_tri_idx = new_tri_indices[j];
                DelunayTriangle& other_tri = triangles[other_tri_idx];

                //Check for shared edge
                int shared_count = 0;
                for (int edge_idx = 1; edge_idx < 3; edge_idx++) {
                    int vi_1 = new_tri.vert_indices[DelunayTriangle::edge_vert_indices[edge_idx][1]];
                    int vi_0 = new_tri.vert_indices[DelunayTriangle::edge_vert_indices[edge_idx][0]];
                    int other_face_idx = other_tri.find_edge(vi_1, vi_0);

                    if (other_face_idx != -1) {
                        //Shared edge
                        new_tri.neighbors[edge_idx] = other_tri_idx;
                        other_tri.neighbors[other_face_idx] = new_tri_idx;
                    }
                }
            }

            new_tri_indices.push_back(new_tri_idx);
        }

    }
}
