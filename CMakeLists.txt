# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.10)
project (CyclopsTetrahedralizer
    VERSION 1.0 
    LANGUAGES CXX
)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()


# Set the C++ standard globally
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_EXECUTABLE "Build executable" ON)
option(BUILD_LIBRARY "Build library" OFF)


file(GLOB_RECURSE sources src/*.cpp src/*.h)

if(BUILD_LIBRARY)
  add_library(cyclopsTetrahedralizerLib STATIC ${sources})
  #Add compile time property
  target_compile_definitions(cyclopsTetrahedralizerLib PRIVATE -DCYCLOPS_TESS_LIB)
endif()

if(BUILD_EXECUTABLE)
  add_executable(cyclopsTetrahedralizer ${sources} "src/delunay_triangulation.h")
endif()



