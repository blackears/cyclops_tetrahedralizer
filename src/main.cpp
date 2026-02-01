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


#include <iostream>

#include <map>
#include <string>

#include "cyclops_tetrahedralizer.h"
#include "obj_file_loader.h"

using namespace std;
using namespace CyclopsTetra3D;

#ifndef CYCLOPS_TESS_LIB

bool parse_command_line(int argc, char** argv, std::map<std::string, std::string>& out_options, std::string& out_source_file)
{
	bool specified_source_file = false;

	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if (arg[0] == '-')
		{
			std::string option_name;
			if (arg[1] == '-') {
				option_name = arg.substr(2);
			}
			else {
				option_name = arg.substr(1);
			}

			std::string option_value;
			if ((i + 1) < argc && argv[i + 1][0] != '-')
			{
				option_value = argv[i + 1];
				++i;
				out_options[option_name] = option_value;
			}

			out_options[option_name] = option_value;
		}
		else
		{
			out_source_file = arg;
			specified_source_file = true;
		}
	}

	if (!specified_source_file)
	{
		cout << "No source file specified." << endl;
		return false;
	}
	return true;
}

void print_help(bool full = false) {
	cout << "Usage" << endl;
	cout << endl;
	cout << "\tcyclopsTetrahedralizer [options] <path to source file>" << endl;
	cout << endl;
	cout << "Specify a source file in the .obj format to read." << endl;

	if (full) {
		cout << endl;
		cout << "\t-o, --out       output .obj file that will be written" << endl;
		cout << "\t-h, --help      help message" << endl;
	}
	else {
		cout << endl;
		cout << "Run 'cyclopsTetrahedralizer --help' for more information." << endl;
	}
}

int main(int argc, char **argv)
{
	std::map<std::string, std::string> options;
	std::string source_file;

	if (!parse_command_line(argc, argv, options, source_file))
	{
		print_help();
		return 1;
	}

	if (options.find("h") != options.end() || options.find("help") != options.end()) {
		print_help(true);
		return 0;
	}

	std::string output_file;
	if (options.find("o") != options.end()) {
		output_file = options["o"];
	}
	else if (options.find("out") != options.end()) {
		output_file = options["out"];
	}
	else {
		int dot_idx = source_file.find_last_of(".");
		if (dot_idx == std::string::npos)
			output_file = source_file + "_tetra.obj";
		else
			output_file = source_file.substr(0, dot_idx) + "_tetra.obj";
	}


	// Load the source file
	WavefrontObjFile loader;
	if (!loader.load_obj_file(source_file))
	{
		cout << "Failed to load source file: " << source_file << endl;
		return 1;
	}

	std::vector<int> face_vertex_indices;
	loader.triangularized_indices(face_vertex_indices);

	//Tetrahedralize
	CyclopsTetrahedralizer tetralizer;
	tetralizer.create_tetrahedrons(loader.get_points(), face_vertex_indices);

	//Export mesh
	//std::vector<int> tri_mesh_vert_indices;
	//tetralizer.get_tetrahedra_as_tri_mesh_indices(tri_mesh_vert_indices);
	//std::vector<int> tri_mesh_face_vert_counts;
	//tri_mesh_face_vert_counts.resize(tri_mesh_vert_indices.size() / 3);
	//std::fill(tri_mesh_face_vert_counts.begin(), tri_mesh_face_vert_counts.end(), 3);
	//WavefrontObjFile result(tetralizer.get_points(), tri_mesh_vert_indices, tri_mesh_face_vert_counts);


	tetralizer.save_obj_file(output_file);

	cout << "Writing tetrahedralization: " << output_file << endl;
	return 0;
}

#endif
