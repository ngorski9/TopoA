cmake_minimum_required(VERSION 3.20)

project(topoCompression LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(VTK
  COMPONENTS
  IOXML
)

add_executable(TopoA src/cli.cpp
                   src/evtk.cpp
                   src/unionfind.cpp
                   src/multitracker.cpp
                   src/unorderedunionfind.cpp)

target_link_libraries(TopoA
  PRIVATE
   ${VTK_LIBRARIES}
   stdc++fs
)

vtk_module_autoinit(
  TARGETS TopoA
  MODULES ${VTK_LIBRARIES}
)

message("it built somehow?")
