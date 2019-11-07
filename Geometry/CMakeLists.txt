cmake_minimum_required(VERSION 3.1)

project(AndersonADMM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
	
	
add_executable(PlanarityOpt
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	AndersonAcceleration.h
	GeometrySolver.h
    Parameters.h
	PlanarityOpt.cpp
	)
	
add_executable(WireMeshOpt
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	AndersonAcceleration.h
	GeometrySolver.h
    Parameters.h
	WireMeshOpt.cpp
	)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message("Clang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP ON CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
	message("AppleClang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP OFF CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message("GNU compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP ON CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message("MSVC compiler found.")
	add_definitions(/DUSE_MSVC)
	add_definitions(/D_USE_MATH_DEFINES)
	add_definitions(/DNOMINMAX)
	set(OPENMP ON CACHE BOOL "OpenMP")
endif()




# Add the OpenMesh folder into include path
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/external")

# OpenMesh library
set(OPENMESH_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/OpenMesh/Core")
file(GLOB_RECURSE OPENMESH_HEADERS "${OPENMESH_SOURCE_DIR}/*.hh" "${OPENMESH_SOURCE_DIR}/*.h")
file(GLOB_RECURSE OPENMESH_SOURCES "${OPENMESH_SOURCE_DIR}/*.cc")
add_library(OpenMesh STATIC ${OPENMESH_HEADERS} ${OPENMESH_SOURCES})
add_definitions(-DOM_STATIC_BUILD)

target_link_libraries(PlanarityOpt OpenMesh)
target_link_libraries(WireMeshOpt OpenMesh)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/Eigen/Dense)
	message("Found Eigen in local folder")
	set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")
	target_include_directories(PlanarityOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(WireMeshOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
else()
	find_package(Eigen3 REQUIRED)
	if(EIGEN3_FOUND)
		target_include_directories(PlanarityOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		target_include_directories(WireMeshOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	endif()
endif()


#OpenMP
set(OPENMP ON CACHE BOOL "OpenMP")
if(OPENMP)
  FIND_PACKAGE(OpenMP QUIET)
  if(OPENMP_FOUND)    
      target_compile_options(PlanarityOpt PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(PlanarityOpt PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(PlanarityOpt "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")  
      target_compile_options(WireMeshOpt PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(WireMeshOpt PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(WireMeshOpt "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      add_definitions(-DEIGEN_DONT_PARALLELIZE)
  else()
      message("OpenMP not found.")
  endif()
endif()

if(UNIX AND NOT APPLE)
find_package(Threads REQUIRED)
target_link_libraries(WireMeshOpt ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(PlanarityOpt ${CMAKE_THREAD_LIBS_INIT})
endif()
