Based on [Bowyer–Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm#)


## Compiling

Get [CMake](https://cmake.org/)


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

### References

Tetrahedralizer plugin for Blender
Matthias Muller
(https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py)

How to build a BVH
jbikker
(https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/)

Bounding Box Ray Intersection
(https://en.wikipedia.org/wiki/Slab_method)