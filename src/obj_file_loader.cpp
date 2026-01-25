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

#include "obj_file_loader.h"

#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <regex>

using namespace CyclopsTetra3D;

bool ObjFileLoader::load_obj_file(const std::string& filename) {
    points.clear();
    face_vertex_indices.clear();
    face_vertex_counts.clear();

    std::ifstream file(filename);
    const std::string whitespace = " \n\r\t\f\v";

    const std::regex command_regex(R"(^(\w+)\s*(.*)$)");
    //const std::regex ws_regex("\\s+");
    const std::regex verts_regex(R"(^\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?))(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?))\s*$)");

    std::string line;
    if (!file.is_open())
        return false;

    while (std::getline(file, line)) {
        if (line.empty())
            continue;

        std::smatch match;
        if (!std::regex_match(line, match, command_regex))
            continue;

        std::string command = match[1];
        std::string args = match[2];

        if (command == "v") {
            std::smatch vert_match;
            if (std::regex_match(args, vert_match, verts_regex)) {
                real x = std::stof(vert_match[1]);
                real y = std::stof(vert_match[4]);
                real z = std::stof(vert_match[7]);
                points.emplace_back(x, y, z);
            }

        }
        else if (command == "f") {
            //std::vector<int> face_inds;
            size_t pos = 0;
            int vertex_count = 0;
            while (pos < args.length()) {
                // Skip leading whitespace
                pos = args.find_first_not_of(whitespace, pos);
                if (pos == std::string::npos)
                    break;

                // Find the end of the current index
                size_t end_pos = args.find_first_of(whitespace, pos);
                std::string index_str = args.substr(pos, end_pos - pos);

                // Extract the vertex index (before any slashes)
                size_t slash_pos = index_str.find('/');
                std::string vertex_index_str = (slash_pos == std::string::npos) ? index_str : index_str.substr(0, slash_pos);

                try {
                    int vertex_index = std::stoi(vertex_index_str);
                    //face_inds.push_back(vertex_index - 1); // OBJ indices are 1-based
                    face_vertex_indices.push_back(vertex_index - 1);
                }
                catch (const std::invalid_argument&) {
                    // Handle invalid index
                    face_vertex_indices.push_back(0);
                }

                vertex_count++;

                if (end_pos == std::string::npos)
                    break;
                pos = end_pos + 1;

            }
            
            face_vertex_counts.push_back(vertex_count);
        }
    }
    
    return true;
}


void ObjFileLoader::save_obj_file(const std::string& filename, const std::vector<Vector3>& points) {
    FILE* file = nullptr;
    errno_t err = fopen_s(&file, filename.c_str(), "w");
    if (file == nullptr) {
        throw std::runtime_error("Failed to open OBJ file for writing: " + filename);
    }

    for (const auto& p : points) {
        fprintf(file, "v %f %f %f\n", p.x, p.y, p.z);
    }

    int face_vertex_index = 0;
    for (int index_count : face_vertex_counts) {
        fprintf(file, "f");
        for (int i = 0; i < index_count; ++i) {
            fprintf(file, " %d", face_vertex_indices[face_vertex_index++] + 1); // OBJ indices are 1-based
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

//void ObjFileLoader::triangularize_face(int face_index_start, int face_vert_count, std::vector<int>& tri_face_indices) {
//    int best_index = 0;
//    real best_radius = INFINITY;
//
//    for (int i = 0; i < face_vert_count; ++i) {
//
//    }
//
//    //if (face_vert_count == 3) {
//    //    tri_face_indices.push_back(face_vertex_indices[face_index_start]);
//    //    tri_face_indices.push_back(face_vertex_indices[face_index_start + 1]);
//    //    tri_face_indices.push_back(face_vertex_indices[face_index_start + 2]);
//    //}
//    //else if (face_vert_count == 4) {
//    //    int fvi0 = face_vertex_indices[face_index_start];
//    //    int fvi1 = face_vertex_indices[face_index_start + 1];
//    //    int fvi2 = face_vertex_indices[face_index_start + 2];
//    //    int fvi3 = face_vertex_indices[face_index_start + 3];
//    //    
//    //    Vector3 p0 = points[fvi0];
//    //    Vector3 p1 = points[fvi1];
//    //    Vector3 p2 = points[fvi2];
//    //    Vector3 p3 = points[fvi3];
//
//    //    Vector3 center_012 = (p0 + p1 + p2) / 3.0;
//    //    Vector3 center_023 = (p0 + p2 + p3) / 3.0;
//    //    Vector3 center_013 = (p0 + p1 + p3) / 3.0;
//    //    Vector3 center_312 = (p1 + p2 + p3) / 3.0;
//    //    
//    //    real split_02 = std::min(
//    //        std::min((p0 - center_023).magnitude_squared(), 
//    //        (p2 - center_023).magnitude_squared(),
//    //        (p3 - center_023).magnitude_squared()),
//
//    //        std::min((p0 - center_012).magnitude_squared(),
//    //        (p1 - center_012).magnitude_squared(),
//    //        (p2 - center_012).magnitude_squared())
//    //    );
//    //    real split_13 = std::min(
//    //        std::min((p0 - center_013).magnitude_squared(),
//    //            (p1 - center_013).magnitude_squared(),
//    //            (p3 - center_013).magnitude_squared()),
//
//    //        std::min((p3 - center_312).magnitude_squared(),
//    //            (p1 - center_312).magnitude_squared(),
//    //            (p2 - center_312).magnitude_squared())
//    //    );
//    //    
//    //    if (split_02 > split_13) {
//
//    //    }
//    //}
//}

void ObjFileLoader::triangularize() {
    std::vector<int> new_face_vertex_indices;
    std::vector<int> new_face_vertex_counts;


    int face_ptr = 0;

    int face_vtx_cursor = 0;
    for (int face_size : face_vertex_counts) {
        if (face_size == 3) {
            new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor]);
            new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 1]);
            new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 2]);
            new_face_vertex_counts.push_back(3);
            face_vtx_cursor += 3;
        }
        else if (face_size > 3) {
            Vector3 center;
            for (int i = 0; i < face_size; ++i) {
                center += points[face_vertex_indices[face_vtx_cursor + i]];
            }
            center /= face_size;

            Vector3 p0 = points[face_vertex_indices[face_vtx_cursor]];
            Vector3 p1 = points[face_vertex_indices[face_vtx_cursor + 1]];
            Vector3 normal = (p0 - center).cross(p1 - center).normalized();

            Plane proj_plane(normal, center);
            Vector3 tangent_0 = normal.cross(Vector3::Y_POS);
            Vector3 tangent_1 = normal.cross(Vector3::X_POS);
            Vector3 tangent = tangent_0.magnitude_squared() > tangent_1.magnitude_squared() ? tangent_0.normalized() : tangent_1.normalized();
            
            Vector3 binormal = tangent.cross(normal);

            std::vector<Vector2> face_points;
            face_points.reserve(face_size);

            for (int i = 0; i < face_size; ++i) {
                Vector3 p = points[face_vertex_indices[face_vtx_cursor + i]];
                Vector2 p_proj = Vector2(tangent.dot(p), binormal.dot(p));
                face_points.push_back(p_proj);
            }

            //std::vector<int> tess_indices = delunay_tessellation(face_points);


        }
        face_ptr += face_size;
    }

    // int face_vtx_cursor = 0;
    // for (int face_idx = 0; face_idx < face_vertex_counts.size(); ++face_idx) {

    //     int vert_count = face_vertex_counts[face_idx];
    //     for (int vert_idx = 2; vert_idx < vert_count; ++vert_idx) {
    //         if ((vert_idx & 0x1) == 0) {
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor]);
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 1]);
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 2]);
    //             new_face_vertex_counts.push_back(3);
    //         }
    //         else {
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor]);
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 2]);
    //             new_face_vertex_indices.push_back(face_vertex_indices[face_vtx_cursor + 1]);
    //             new_face_vertex_counts.push_back(3);
    //         }
    //         face_vtx_cursor += 1;
    //     }
    // }

    
    face_vertex_indices = std::move(new_face_vertex_indices);
    face_vertex_counts = std::move(new_face_vertex_counts);
}
