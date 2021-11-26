cmake_minimum_required(VERSION 3.20)
project(openMeshPolyScope)

set(CMAKE_CXX_STANDARD 14)

include_directories(deps)
set(open_mesh_sources
        deps/OpenMesh
        deps/OpenMesh/Core
        deps/OpenMesh/Core/Geometry
        deps/OpenMesh/Core/IO
        deps/OpenMesh/Core/IO/exporter
        deps/OpenMesh/Core/IO/importer
        deps/OpenMesh/Core/IO/reader
        deps/OpenMesh/Core/IO/writer
        deps/OpenMesh/Core/Mesh
        deps/OpenMesh/Core/Mesh/gen
        deps/OpenMesh/Core/System
        deps/OpenMesh/Core/Templates
        deps/OpenMesh/Core/Utils
        )

# Add headers recursively to headers
file(GLOB_RECURSE om_headers "*.h*" ${open_mesh_sources})

# Add sources recursively to sources
file(GLOB_RECURSE om_sources "*.cc" ${open_mesh_sources})

# We need a static build of Open mesh
add_definitions(-DOM_STATIC_BUILD)

# Needed for windows compilation
if (WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif ()

# Polyscope
add_subdirectory(deps/polyscope)

add_executable(openMeshPolyScope main.cpp
        ${om_headers}
        ${om_sources})

target_link_libraries(openMeshPolyScope polyscope)
