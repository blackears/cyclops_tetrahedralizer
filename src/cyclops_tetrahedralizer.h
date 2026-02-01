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


#ifndef CYCLOPS_TETRAHEDRALIZER_H
#define CYCLOPS_TETRAHEDRALIZER_H

#include <vector>
#include <array>

#include "math.h"

namespace CyclopsTetra3D {



//Face winding - face normals points outward
class Tetrahedron {
    //Vertex ordering per face
    //Tetrahedron faces wind ccw when viewed from outside
    static constexpr std::array<std::array<int, 3>, 4> face_vert_indices = { {
        {0, 1, 2}, {1, 0, 3}, {2, 3, 0}, {3, 2, 1}
        } };

private:
    std::array<int, 4> vert_indices;
    std::array<int, 4> neighbors;

    Vector3 circumcenter;
    real circumcircle_radius_squared;
    Vector3 center;

    std::array<Plane, 4> face_planes;

    bool valid;
public:

    static Tetrahedron create_from_points(int v0_idx, int v1_idx, int v2_idx, int v3_idx, const std::vector<Vector3>& points);

    const std::array<int, 4>& get_vert_indices() const { return vert_indices; }

    bool has_neighbor(int neighbor_idx) const {
        for (int i = 0; i < 4; i++) {
            if (neighbors[i] == neighbor_idx) {
                return true;
            }
        }
        return false;
    }
    
    bool point_in_circumsphere(const Vector3& p, const std::vector<Vector3>& points) const {
        return (circumcenter - p).magnitude_squared() < circumcircle_radius_squared;
    }
    bool contains_point(const Vector3& p, const std::vector<Vector3>& points) const;
    int find_adjacent_tetrahedron(const Vector3& dir, const std::vector<Vector3>& points) const;
    real quality(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) const;

    //vertex indices must wind face outward
    int find_face(int v0_idx, int v1_idx, int v2_idx) const {
        for (int i = 0; i < 4; i++) {
            if ((vert_indices[face_vert_indices[i][0]] == v0_idx &&
                 vert_indices[face_vert_indices[i][1]] == v1_idx &&
                 vert_indices[face_vert_indices[i][2]] == v2_idx) ||
                (vert_indices[face_vert_indices[i][0]] == v1_idx &&
                 vert_indices[face_vert_indices[i][1]] == v2_idx &&
                 vert_indices[face_vert_indices[i][2]] == v0_idx) ||
                (vert_indices[face_vert_indices[i][0]] == v2_idx &&
                 vert_indices[face_vert_indices[i][1]] == v0_idx &&
                 vert_indices[face_vert_indices[i][2]] == v1_idx)) {
                return i;
            }
        }
        return -1;
    }

    friend class CyclopsTetrahedralizer;
};

class CyclopsTetrahedralizer {
    std::vector<Vector3> tess_points;
    std::vector<Tetrahedron> tetrahedra;

private:
    void create_tetrahedrons_iter(std::vector<Tetrahedron>& tetrahedrons, const std::vector<Vector3>& points);

public:
    const std::vector<Vector3>& get_points() const { return tess_points; }
    std::vector<Vector3>& get_points() { return tess_points; }
    const std::vector<Tetrahedron>& get_tetrahedra() const { return tetrahedra; }

    int num_valid_tetrahedra() const {
        int count = 0;
        for (auto& tet : tetrahedra) {
            if (tet.valid)
                count++;
        }
        return count;
    }

    const void get_tetrahedra_as_tri_mesh_indices(std::vector<int>& out_tri_indices) const { 
        out_tri_indices.resize(num_valid_tetrahedra() * 12);
        int offset = 0;
        for (const Tetrahedron& tet : tetrahedra) {
            if (!tet.valid)
                continue;

            for (int i = 0; i < 4; ++i) {
                out_tri_indices[offset++] = tet.vert_indices[Tetrahedron::face_vert_indices[i][0]];
                out_tri_indices[offset++] = tet.vert_indices[Tetrahedron::face_vert_indices[i][2]];
                out_tri_indices[offset++] = tet.vert_indices[Tetrahedron::face_vert_indices[i][3]];
            }
        }
    }

    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    //@param resolution spacing for extra interior points
    void create_tetrahedrons(const std::vector<Vector3>& points, const std::vector<int>& indices, float resolution = 0);

    void get_mesh(std::vector<Vector3>& out_points, std::vector<int>& out_indices);

    void save_obj_file(const std::string& filename) const;

    };

}

#endif //CYCLOPS_TETRAHEDRALIZER_H