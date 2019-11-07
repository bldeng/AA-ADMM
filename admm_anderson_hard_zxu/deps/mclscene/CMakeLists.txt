# Copyright (c) 2017 University of Minnesota
# 
# MCLSCENE Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
#    of conditions and the following disclaimer in the documentation and/or other materials
#    provided with the distribution.
# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# By Matt Overby (http://www.mattoverby.net)

cmake_minimum_required(VERSION 3.1)
project(mclscene C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_BUILD_TYPE Release)
add_definitions( -DMCLSCENE_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" )

# MCL options
option(MCL_BUILD_TESTS "Build MCL tests" ON)
option(MCL_BUILD_EXAMPLES "Build MCL examples" ON)
option(MCL_USE_GLEW "Build with GLEW" ON)

# Compiler options
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long")
endif()

############################################################
#
#	Libraries
#
############################################################

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Eigen3 (header only)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# mclscene (header only), submodule in deps/
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# glfw3
find_package(GLFW 3 REQUIRED)
include_directories(SYSTEM ${GLFW_INCLUDE_DIR})

# glew
if(MCL_USE_GLEW)
add_definitions(-DMCL_USE_GLEW)
find_package(GLEW REQUIRED)
include_directories(SYSTEM ${GLEW_INCLUDE_DIR})
else()
	set(GLEW_LIBRARIES "")
endif()

# opengl
find_package(OpenGL 3.2 REQUIRED)
include_directories(SYSTEM ${OPENGL_INCLUDE_DIR})

# var to hold the libs
set(RENDER_LIBS
	${OPENGL_LIBRARIES}
	${GLEW_LIBRARIES}
	${GLFW_LIBRARIES}
)

# tetgen, creates library "tet"
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/tetgen)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/src/tetgen)

#
# ILMBase, required for Alembic
#
FIND_PACKAGE(IlmBase)
option(MCL_BUILD_ABC_EXAMPLES "Compile Alembic examples" on)
if(ILMBASE_FOUND)
	include_directories(SYSTEM ${ILMBASE_INCLUDE_DIRECTORY} )
endif()

#
# Alembic: see http://www.alembic.io/
#
find_package(Alembic)
if(ALEMBIC_FOUND AND ILMBASE_FOUND)
	include_directories(SYSTEM ${ALEMBIC_INCLUDE_PATH}/../ )
else()
	# Turn off build of alembic examples if ILMBase/Alembic libs not found
	set(MCL_BUILD_ABC_EXAMPLES off)
endif()


#
# I know I should just use add_subdirectory and keep these
# in a seperate CMakeLists.txt. However, it's more convenient to keep
# them all in one while the project is still early in development.
#

# Build tests
if(MCL_BUILD_TESTS)

	enable_testing()

	add_executable(test_graphcolor src/tests/test_graphcolor.cpp)
	add_test(test_graphcolor test_graphcolor)

	add_executable(test_embeddedmesh src/tests/test_embeddedmesh.cpp)
	add_test(test_embeddedmesh test_embeddedmesh)

	add_executable(test_tetgen src/tests/test_tetgen.cpp)
	target_link_libraries(test_tetgen tet)
	add_test(test_tetgen test_tetgen)

	add_executable(test_varmanager src/tests/test_varmanager.cpp)
	add_test(test_varmanager test_varmanager)

	add_executable(test_timer src/tests/test_timer.cpp)
	add_test(test_timer test_timer)

endif(MCL_BUILD_TESTS)

# Build examples
if(MCL_BUILD_EXAMPLES)

	add_executable(renderWindow src/examples/renderWindow.cpp)
	target_link_libraries(renderWindow ${RENDER_LIBS})

	add_executable(makeTetmesh src/examples/makeTetmesh.cpp)
	target_link_libraries(makeTetmesh tet ${RENDER_LIBS})

	add_executable(tetmeshViewer src/examples/tetmeshViewer.cpp)
	target_link_libraries(tetmeshViewer ${RENDER_LIBS})

	if(MCL_BUILD_ABC_EXAMPLES)
		add_executable(alembicExport src/examples/alembicExport.cpp)
		target_link_libraries(alembicExport ${RENDER_LIBS} ${ILMBASE_LIBS} ${ALEMBIC_LIB} )

		add_executable(arcsimToAlembic src/examples/arcsimToAlembic.cpp)
		target_link_libraries(arcsimToAlembic ${RENDER_LIBS} ${ILMBASE_LIBS} ${ALEMBIC_LIB} )
	endif()

endif(MCL_BUILD_EXAMPLES)


