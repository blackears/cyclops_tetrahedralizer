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


#ifndef DELUNAY_TRIANGULATION_H
#define DELUNAY_TRIANGULATION_H

#include  <vector>
#include "math.h"

namespace CyclopsTetra3D {

//Face winding - face normals points outward
struct DelunayTriangle {
    static constexpr int edge_vert_indices[3][2] = { {0, 1}, {1, 2}, {2, 0} };

    int vert_indices[3];
    int neighbors[3];

    Vector2 circumcenter;
    real circumcircle_radius_squared;
    Vector2 center;

    bool valid;

    void create_from_points(int v0_idx, int v1_idx, int v2_idx, const std::vector<Vector2>& points) {
        vert_indices[0] = v0_idx;
        vert_indices[1] = v1_idx;
        vert_indices[2] = v2_idx;

        neighbors[0] = -1;
        neighbors[1] = -1;
        neighbors[2] = -1;

        circumcenter = Math::triangle_circumcenter(points[v0_idx], points[v1_idx], points[v2_idx]);
        circumcircle_radius_squared = (circumcenter - points[vert_indices[0]]).magnitude_squared();
        center = (points[v0_idx] + points[v1_idx] + points[v2_idx]) / 3.0;

        valid = true;
    }

    bool has_neighbor(int neighbor_idx) const {
        for (int i = 0; i < 3; i++) {
            if (neighbors[i] == neighbor_idx) {
                return true;
            }
        }
        return false;
    }
    bool point_in_circumcircle(const Vector2& p, const std::vector<Vector2>& points) const {
        return (circumcenter - p).magnitude_squared() < circumcircle_radius_squared;
    }
    bool contains_point(const Vector2& p, const std::vector<Vector2>& points) const {
        return Math::triangle_contains_point(p, points[vert_indices[0]], points[vert_indices[1]], points[vert_indices[2]]);
    }

    int find_adjacent_triangle(const Vector2& dir, const std::vector<Vector2>& points) const {
        Vector2 p0 = points[vert_indices[0]];
        Vector2 p1 = points[vert_indices[1]];
        Vector2 p2 = points[vert_indices[2]];

        Vector2 hit_0 = Math::intersect_lines(center, dir, p0, p1 - p0);
        Vector2 hit_1 = Math::intersect_lines(center, dir, p0, p2 - p0);
        Vector2 hit_2 = Math::intersect_lines(center, dir, p1, p2 - p1);

        real dot_0 = hit_0.dot(dir);
        real dot_1 = hit_1.dot(dir);
        real dot_2 = hit_2.dot(dir);

        //Smallest positive dot product
        if (dot_0 >= 0.0 && (dot_1 < 0.0 || dot_1 > dot_0) && (dot_2 < 0.0 || dot_2 > dot_0))
            return neighbors[0];
        if (dot_1 >= 0.0 && (dot_2 < 0.0 || dot_2 > dot_0))
            return neighbors[1];
        return neighbors[2];
    }

    int find_edge(int v0_idx, int v1_idx) const {
        if (v0_idx == vert_indices[0] && v1_idx == vert_indices[1])
            return 0;
        if (v0_idx == vert_indices[1] && v1_idx == vert_indices[2])
            return 1;
        if (v0_idx == vert_indices[2] && v1_idx == vert_indices[0])
            return 2;
        return -1;
    }

};


class DelunayTriangulator {
    std::vector<DelunayTriangle> triangles;

private:
    void create_triangles_iter(const std::vector<Vector2>& points);

public:
    void create_triangles(const std::vector<Vector2>& points, const std::vector<int>& indices, float resolution = 0);

    const std::vector<DelunayTriangle>& get_triangles() const { return triangles; }
};

}

#endif