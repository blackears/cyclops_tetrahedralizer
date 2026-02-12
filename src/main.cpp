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

#include <algorithm>
#include <string>

#include "cyclops_tetrahedralizer.h"
#include "wavefront_obj_file.h"

using namespace std;
using namespace CyclopsTetra3D;

#ifndef CYCLOPS_TESS_LIB

bool show_help = false;
bool export_edges = false;
real subdivisions = 0.0;
std::string output_file;

void process_option(std::string option, int argc, char** argv, int& option_ptr) {
	if (option == "h" || option == "help") {
		show_help = true;
		return;
	}

	if (option == "e" || option == "edges") {
		export_edges = true;
		return;
	}

	if (option == "o" || option == "out") {
		output_file = argv[option_ptr++];
		return;
	}

	if (option == "s" || option == "subdiv") {
		char* last_char;
		double result = std::strtod(argv[option_ptr++], &last_char);
		if (*last_char == 0) {
			subdivisions = std::max(result, 0.0);
		}
		return;
	}
}


bool parse_command_line(int argc, char** argv, std::string& out_source_file)
{
	bool specified_source_file = false;

	for (int i = 1; i < argc;)
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
			i++;

			process_option(option_name, argc, argv, i);
		}
		else
		{
			out_source_file = arg;
			specified_source_file = true;
			i++;
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
		cout << "\t-h, --help                  help message" << endl;
		cout << "\t-o, --out <filename>        output .obj file that will be written.  If not" << endl;
		cout << "\t                            specified, a file will be written where the suffix" << endl;
		cout << "\t                            \"_tetra\" is appended to the input file name." << endl;
		cout << "\t-e, --edges                 export edges instead of faces" << endl;
		cout << "\t-s, --subdiv <number>       if greater than 0, applies a cube grid to" << endl;
		cout << "\t                            mesh, with cube size is the length of the max" << endl;
		cout << "\t                            side legth of the bounding box divided into" << endl;
		cout << "\t                            this many segments" << endl;
	}
	else {
		cout << endl;
		cout << "Run 'cyclopsTetrahedralizer --help' for more information." << endl;
	}
}

int main(int argc, char **argv)
{
	//std::map<std::string, std::string> options;
	std::string source_file;

	if (!parse_command_line(argc, argv, source_file))
	{
		print_help();
		return 1;
	}

	if (show_help) {
		print_help(true);
		return 0;
	}

	if (output_file.empty()) {
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
	tetralizer.create_tetrahedrons(loader.get_points(), face_vertex_indices, subdivisions);

	cout << "tessellation done" << std::endl;

	if (export_edges)
		tetralizer.save_file_line_segments_obj(output_file);
	else {
		tetralizer.save_file_obj(output_file);
	}

	cout << "Writing tetrahedralization: " << output_file << endl;
	return 0;
}

#endif
