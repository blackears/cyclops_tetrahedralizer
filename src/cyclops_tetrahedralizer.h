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

#include  <vector>

#include "math.h"

namespace CyclopsTetra3D {



//Face winding - face normals points outward
struct Tetrahedron {
    //Vertex ordering per face
    static constexpr int face_vert_indices[4][3] = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {0, 2, 3}};
    static constexpr int face_missing_vert_index[4] = {3, 2, 0, 1};

    int vert_indices[4];
    int neighbors[4];

    Vector3 circumcenter;
    real circumcircle_radius_squared;
    Vector3 center;

    Plane face_planes[4];

    bool valid;

    static Tetrahedron create_from_points(int v0_idx, int v1_idx, int v2_idx, int v3_idx, const std::vector<Vector3>& points);

    bool has_neighbor(int neighbor_idx) const {
        for (int i = 0; i < 4; i++) {
            if (neighbors[i] == neighbor_idx) {
                return true;
            }
        }
        return false;
    }
    
    //static Vector3 calc_circumcenter(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) {
    //    //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter
    //    
    //    //From Matthias Muller
    //    //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68
    //    Vector3 b = p1 - p0;
    //    Vector3 c = p2 - p0;
    //    Vector3 d = p3 - p0;

    //    real det = 2.0 * (b.x * (c.y * d.z - c.z * d.y)
    //        - b.y * (c.x * d.z - c.z * d.x)
    //        + b.z * (c.x * d.y - c.y * d.x));

    //    if (det == 0.0) {
    //        return p0;
    //    }
    //    else {
    //        Vector3 v = c.cross(d) * b.dot(b) + d.cross(b) * c.dot(c) + b.cross(c) * d.dot(d);
    //        v /= det;
    //        return p0 + v;
    //    }
    //}

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
};

class CyclopsTetrahedralizer {
    Vector3* point_list;
    std::vector<Tetrahedron> tetrahedrons;

private:
    void create_tetrahedrons_iter(std::vector<Tetrahedron>& tetrahedrons, const std::vector<Vector3>& points);

public:
    //@param points of triangles
    //@param indices of triangles (3 per triangle)
    //@param resolution spacing for extra interior points
    void create_tetrahedrons(const std::vector<Vector3>& points, const std::vector<int>& indices, float resolution = 0);

};

}

#endif //CYCLOPS_TETRAHEDRALIZER_H