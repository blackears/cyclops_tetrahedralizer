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
			std::string option_name = arg.substr(1);
			std::string option_value;
			if ((i + 1) < argc)
			{
				option_value = argv[i + 1];
				++i;
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

int main(int argc, char **argv)
{
	std::map<std::string, std::string> options;
	std::string source_file;

	if (!parse_command_line(argc, argv, options, source_file))
	{
		return 1;
	}

	// Load the source file
	ObjFileLoader<float> loader;
	loader.load_obj_file(source_file);
	// if (!loader.load_obj_file(source_file))
	// {
	// 	cout << "Failed to load source file: " << source_file << endl;
	// 	return 1;
	// }

	cout << "Hello CMake3." << endl;
	return 0;
}

#endif
