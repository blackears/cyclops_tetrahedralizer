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
    face_indices.clear();
    face_index_counts.clear();

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
                    face_indices.push_back(vertex_index);
                }
                catch (const std::invalid_argument&) {
                    // Handle invalid index
                    face_indices.push_back(0);
                }

                vertex_count++;

                if (end_pos == std::string::npos)
                    break;
                pos = end_pos + 1;

            }
            
            face_index_counts.push_back(vertex_count);
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
    for (int index_count : face_index_counts) {
        fprintf(file, "f");
        for (int i = 0; i < index_count; ++i) {
            fprintf(file, " %d", face_indices[face_vertex_index++]); // OBJ indices are 1-based
        }
        fprintf(file, "\n");
    }

    fclose(file);
}
