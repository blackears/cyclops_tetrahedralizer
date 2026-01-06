// Cyclops Tessellate.cpp : Defines the entry point for the application.
//

#include <iostream>

using namespace std;

#ifdef CYCLOPS_TESS_LIB
void tessellate_tetrahedra(float* points) {

#else

int main()
{
	cout << "Hello CMake2." << endl;
	return 0;
}


#endif