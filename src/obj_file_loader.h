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

/*
 * In .obj files, clockwise winding triangles point outward.
 */
 
#ifndef CYCLOPS_OBJ_FILE_LOADER_H
#define CYCLOPS_OBJ_FILE_LOADER_H

#include <vector>
#include <string>

#include "math.h"

namespace CyclopsTetra3D {

class ObjFileLoader
{
    std::vector<Vector3> points;
    std::vector<int> face_vertex_indices;
    std::vector<int> face_vertex_counts;

    //Normal points in direction of CCW winding.
    Vector3 calc_face_normal(int face_index_start, int face_vert_count) {
        switch (face_vert_count) {
        case 4: {
            Vector3 p0 = points[face_vertex_indices[face_index_start]];
            Vector3 p1 = points[face_vertex_indices[face_index_start + 1]];
            Vector3 p2 = points[face_vertex_indices[face_index_start + 2]];
            Vector3 p3 = points[face_vertex_indices[face_index_start + 3]];

            return (p2 - p0).cross(p3 - p1).normalized();
        }
        case 3: {
            Vector3 p0 = points[face_vertex_indices[face_index_start]];
            Vector3 p1 = points[face_vertex_indices[face_index_start + 1]];
            Vector3 p2 = points[face_vertex_indices[face_index_start + 2]];

            return (p1 - p0).cross(p2 - p1).normalized();
        }
        default:
            //Newell's method
            //https://people.eecs.berkeley.edu/~ug/slide/pipeline/assignments/backfacecull.shtml

            //n.x is area of polygon projected onto YZ plane, etc.
            Vector3 normal_sum;
            for (int i = 0; i < face_vert_count; ++i) {
                int j = (i < face_vert_count - 1) ? i + 1 : i - face_vert_count + 1;
                Vector3 p0 = points[face_vertex_indices[face_index_start + i]];
                Vector3 p1 = points[face_vertex_indices[face_index_start + j]];

                normal_sum.x += (p1.z + p0.z) * (p1.y - p0.y);
                normal_sum.y += (p1.x + p0.x) * (p1.z - p0.z);
                normal_sum.z += (p1.y + p0.y) * (p1.x - p0.x);
            }
            return normal_sum.normalized();
        }

    }

public:
    ObjFileLoader() = default;
    ~ObjFileLoader() = default;

    bool load_obj_file(const std::string& filename);
    void save_obj_file(const std::string& filename, const std::vector<Vector3>& points);

    void triangularize();
};


}

#endif //CYCLOPS_OBJ_FILE_LOADER_H