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

using namespace CyclopsTetra3D;

bool ObjFileLoader::load_obj_file(const std::string& filename) {
    points.clear();
    face_indices.clear();

    //std::filesystem::path cp = std::filesystem::current_path();

    FILE* file = nullptr;
    errno_t err = fopen_s(&file, filename.c_str(), "r");
    if (file == nullptr) {
        throw std::runtime_error("Failed to open OBJ file: " + filename);
    }
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == 'v' && line[1] == ' ') {
            real x, y, z;
            if (sscanf_s(line + 2, "%f %f %f", &x, &y, &z) == 3) {
                points.emplace_back(x, y, z);
            }
        }

        if (line[0] == 'f' && line[1] == ' ') {
            std::vector<int> face_indices;
            char* ptr = line + 2;
            while (*ptr) {
                int index;
                if (sscanf_s(ptr, "%d", &index) == 1) {
                    face_indices.push_back(index - 1); // OBJ indices are 1-based
                    while (*ptr && *ptr != ' ') ++ptr;
                }
                else {
                    break;
                }
                while (*ptr == ' ') ++ptr;
            }
        }
    }

    fclose(file);
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

    for (const auto& face : face_indices) {
        fprintf(file, "f");
        for (const auto& index : face) {
            fprintf(file, " %d", index + 1); // OBJ indices are 1-based
        }
        fprintf(file, "\n");
    }

    fclose(file);
}
