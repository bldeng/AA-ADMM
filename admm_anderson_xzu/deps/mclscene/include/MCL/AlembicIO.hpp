// Copyright (c) 2017 University of Minnesota
// 
// MCLSCENE Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list
//    of conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE 
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// By Matt Overby (http://www.mattoverby.net)

#ifndef MCL_ALEMBICIO_H
#define MCL_ALEMBICIO_H

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

namespace mcl {

//
//	A very simplified wrapper for writing out to Alembic (.abc) files.
//	Designed for very shallow scenes.
//	Needs a better design, since a forced exit (no call to destructor)
//	will cause problems with the alembic output file.
//
class AlembicExporter {
public:
	typedef std::shared_ptr<Alembic::AbcGeom::OObject> ObjPtr;
	typedef std::shared_ptr<Alembic::AbcGeom::OPolyMesh> MeshPtr;

	AlembicExporter(){}

	AlembicExporter( int framerate, const std::string &filename );

	// Adds an object to the scene, returns its handle (index)
	inline size_t add_object( std::string name = "" );

	// Sets the data at this frame.
	// Assumes CCW winding order as input. Option to output in CW or CCW is provided.
	inline void add_frame( size_t handle, const float *verts, int n_verts,
		const int *inds, int n_prims, bool output_CCW=false );

private:
	std::shared_ptr<Alembic::AbcGeom::OArchive> archive;
	Alembic::AbcCoreAbstract::chrono_t dt;
	std::vector<ObjPtr> objects;
	std::vector<MeshPtr> meshes;

}; // class AlembicExporter


//
//	Implementation
//


AlembicExporter::AlembicExporter( int framerate, const std::string &filename ) :
	dt(1.f/float(framerate)) {
	using namespace Alembic::AbcGeom;

	archive = std::make_shared<OArchive>( Alembic::AbcCoreOgawa::WriteArchive(), filename );
	archive->addTimeSampling(  TimeSampling(dt,0.0) );

} // end constructor


inline size_t AlembicExporter::add_object( std::string name ){
	using namespace Alembic::AbcGeom;

	if( name.length()==0 ){ name = "object" + std::to_string( objects.size() ); }
	size_t handle = objects.size();

	// Create an object and an associated mesh
	objects.emplace_back(
		std::make_shared<OObject>( Alembic::AbcGeom::OObject( *archive.get(), kTop ), name.c_str() )
	);
	meshes.emplace_back(
		std::make_shared<OPolyMesh>( *objects.back(), "mesh", archive->getTimeSampling(1) )
	);

	return handle;
} // end add object


inline void AlembicExporter::add_frame( size_t handle, const float *verts, int n_verts,
	const int *inds, int n_prims, bool output_CCW ){
	using namespace Alembic::AbcGeom;

	if( handle >= meshes.size() ){
		throw std::runtime_error("**AlembicExporter Error: Handle not found. Call add_object first");
	} 
	MeshPtr meshptr = meshes[handle];
	OPolyMeshSchema &mesh = meshptr->getSchema();

	// "Counts" refers to the number of vertices in each face (tris=3)
	std::vector<int> counts( n_prims, 3 );

	// Need to make a copy of indices to flip winding order.
	std::vector<int> indices(n_prims*3);
	for( int i=0; i<n_prims; ++i ){
		indices[i*3] = inds[i*3];
		if( !output_CCW ){
			indices[i*3+2] = inds[i*3+1];
			indices[i*3+1] = inds[i*3+2];
		} else {
			indices[i*3+1] = inds[i*3+1];
			indices[i*3+2] = inds[i*3+2];
		}
	}

	// Because we don't know the type 
	OPolyMeshSchema::Sample sample(
		V3fArraySample((const V3f*)verts, n_verts),
		Int32ArraySample(indices),
		Int32ArraySample(counts)
	);

	mesh.set(sample);

} // end set frame


} // namespace mcl

#endif
