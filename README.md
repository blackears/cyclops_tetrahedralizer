# Cyclops Tetrahedralizer

This is a C++ library meant for finding tetrahedralizations of input meshes.  This is useful for programs like physics simulations that need to break a complex mesh down into a collection of simplexes.  

It also can be used as a command line utitily to read in Alias/Wavefront .obj files and output a file that contains the tetrahedralization of the mesh.

The quality of the mesh is functional, but not the best at the moment.

The motivation behind writing a new library when there are already several available is twofold.  First, Cyclops Tetrahedralizer is available under the MIT license, as other libraries which publish under GPL family licenses create barriers to incorporating the code into non-GPL projects.  Second, this library is meant to be standalone, as other libraries which require dependencies on huge math libraries make them prohibitive to integrate into other projects.

## Usage

The command line tool can be used to read an input Alias/Wavefront .obj file and output a second file containing the tetrahedralized mesh.

	cyclopsTetrahedralizer [options] <path to source file>

Options:

	-h, --help
		help message
	-o, --out <filename>
		output .obj file that will be written.  If not specified, a file will be written where the suffix "_tetra" is appended to the input file name.
	-e, --edges
		export edges instead of faces
	-s, --subdiv <number>
		if greater than 0, applies a cube grid to mesh, with cube size is the length of the max side legth of the bounding box divided into this many segments

Example usage:

	cyclopsTetrahedralizer -o cylinder_tet.obj -s 5 cylinder.obj

## Compiling

### Requirements

[CMake](https://cmake.org/) is required to build this.

### Building the project

From the download directory

	cd <cyclops tessellator project directory>
	cmake -B build -S .

This will build a solution inside of the `./build` directory.  To create the executable, either open the `build/CyclopsTetrahedralizer.sln` in Visual Studio and select Build All.  Alternately, on the command line you can issue:

	msbuild .\build\CyclopsTetrahedralizer.sln

Either option will compile the executable and place it in `Debug/cyclopsTetrahedralizer.exe`.


### On Windows:

Open CyclopsTetrahedralizer.sln project in Visual Studio and build all.
`CyclopsTetrahedralizer.exe` will appear under /build/Debug


### On Mac:
```
make
```

## References

#### Tetrahedralizer plugin for Blender

Matthias Muller

(https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py)


#### How to build a BVH
jbikker

(https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/)

#### Bounding Box Ray Intersection
(https://en.wikipedia.org/wiki/Slab_method)

#### Bowyer–Watson algorithm
[Bowyer–Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm#)
