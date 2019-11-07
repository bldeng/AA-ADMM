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

#include "MCL/MeshIO.hpp"
#include "MCL/AlembicIO.hpp"
#include "MCL/ArgParser.hpp"
#include "MCL/RenderWindow.hpp"
#include "../pugixml/pugixml.hpp"

using namespace mcl;

void help();
bool file_exists( const std::string &filename );
int num_meshes( const std::string &dir );
int num_obs( const std::string &dir );
std::string pad_leading_zeros( int num_digits, int num );
std::string get_name( const std::string &dir );
bool test_conf( const std::string &dir );


std::string dir;
std::string abc_name;
AlembicExporter exporter;
std::vector< TriangleMesh::Ptr > obstacles;
std::vector< Eigen::AlignedBox<float,3> > obs_aabb;
std::vector<int> obs_handles;
std::vector< TriangleMesh::Ptr > meshes;
std::vector<int> deform_handles;
std::vector< RenderMesh::Ptr > render_meshes;

inline std::string make_screenshot_fn( std::string outdir, int &counter ){
	std::stringstream ss;
	ss << outdir << "/" << std::setfill('0') << std::setw(6) << counter << ".png";
	counter++;
	return ss.str();
}

bool load_frame( int frame_num, bool &stop ){

	int n_meshes = meshes.size();
	int n_obs = obstacles.size();

	// Winding order matters for whatever is reading in the abc file.
	bool ccw_output = false;

	std::string framestr = pad_leading_zeros(6,frame_num);

	// Load all deforming meshes
	for( int m=0; m<n_meshes; ++m ){
		std::string meshstr = pad_leading_zeros(2,m);
		std::string objfile = dir + framestr + '_' + meshstr + ".obj";

		// If the file doesn't exists, we have reached the end of the sim
		if( !file_exists(objfile) ){ 
			stop = true;
			break;
		}

		// Load the deforming mesh
		TriangleMesh *mesh = meshes[m].get();
		if( !meshio::load_obj( mesh, objfile, false, false, false ) ){
			std::cerr << "\n**arcsimeToAlembic Error: Failed to load " << objfile << "\n" << std::endl;
			return false;
		}
		mesh->need_normals();

		// Add it to the exporter
		exporter.add_frame( deform_handles[m], &mesh->vertices[0][0], mesh->vertices.size(),
			&mesh->faces[0][0], mesh->faces.size(), ccw_output );

	} // end loop meshes

	// Load all obstacle transforms
	for( int o=0; o<n_obs && !stop; ++o ){
		std::string meshstr = pad_leading_zeros(2,o);
		std::string xformfile = dir + framestr + "obs" + meshstr + ".txt";
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(xformfile.c_str());
		if( !result ){
			std::cerr << "\n**arcsimeToAlembic Error: Unable to load " << xformfile << std::endl;
			return false;
		}

		TriangleMesh *mesh = obstacles[o].get();
		Eigen::AlignedBox<float,3> aabb = obs_aabb[o];
		Vec3f obs_center = aabb.center();

		pugi::xml_node::iterator node_iter = doc.first_child();
		for( ; node_iter != doc.end(); node_iter++ ){
			pugi::xml_node curr_node = *node_iter;
			std::string name = curr_node.name();
			XForm<float> xf; xf.setIdentity();
			if( name == "rotate" ){
				float angle = curr_node.attribute("angle").as_float() * (180.f/M_PI);
				float x = curr_node.attribute("x").as_float();
				float y = curr_node.attribute("y").as_float();
				float z = curr_node.attribute("z").as_float();

				// Translate obs to origin before rotation
				XForm<float> xf0 = xform::make_trans<float>(-obs_center);
				XForm<float> xf1 = xform::make_rot<float>( angle, Vec3f(x,y,z) );
				XForm<float> xf2 = xform::make_trans<float>(obs_center);
				xf = xf2 * xf1 * xf0;
			}
			else if( name == "scale" ){
				float s = curr_node.attribute("value").as_float();
				xf = xform::make_scale<float>(s,s,s);
			}
			else if( name == "translate" ){
				float x = curr_node.attribute("x").as_float();
				float y = curr_node.attribute("y").as_float();
				float z = curr_node.attribute("z").as_float();
				xf = xform::make_trans<float>(x,y,z);
			}
			mesh->apply_xform(xf);
		} // end load xform

		// Add it to the exporter
		exporter.add_frame( obs_handles[o], &mesh->vertices[0][0], mesh->vertices.size(),
			&mesh->faces[0][0], mesh->faces.size(), ccw_output );

	} // end loop obstacles

	return true;
}

bool framedump_running;
class DumpController : public Controller {
public:
	DumpController(){
		framedump_running = false;
		std::cout << "\nMove the camera to the location you want it and press SPACEBAR to start\n" << std::endl;
	}
	void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
		Controller::key_callback(window,key,scancode,action,mods); // call base class cb if you want
		if( action != GLFW_PRESS ){ return; }
		switch ( key ){
			case GLFW_KEY_SPACE:{
				framedump_running = !framedump_running;
			} break;
		}
	}
};

int main(int argc, char *argv[]){

	ArgParser parser(argc,argv);
	if(	parser.exists("-h") || parser.exists("-help") ||
		parser.exists("--h") || parser.exists("--help") ||
		!parser.exists("-i") || !parser.exists("-dt") ){
		help();
		return EXIT_SUCCESS;
	}

	dir = parser.get<std::string>("-i");
	if( dir.back()!='/' ){ dir += '/'; }
	unsigned int skipframe = 1;
	parser.get<unsigned int>( "-skip", &skipframe );
	double dt = parser.get<double>("-dt");
	if( skipframe > 2 ){ dt += 1.0 / double(skipframe); }
	bool save_frames = false;
	std::string frames_dir = "";
	parser.get<std::string>( "-dump", &frames_dir );
	if( frames_dir.length() > 0 ){ save_frames = true; }

	// Attempt to open conf.json to make sure the diretory is all dandy.
	if( !test_conf(dir) ){ return EXIT_FAILURE; }

	// Create the ABC file
	abc_name = get_name( dir );
	abc_name = abc_name + ".abc";
	parser.get<std::string>("-o", &abc_name);
	exporter = AlembicExporter( int(1.0/dt), abc_name );

	// Get the number of meshes and obs
	int n_meshes = num_meshes( dir );
	if( n_meshes <= 0 ){ return EXIT_FAILURE; }
	int n_obs = num_obs( dir );
	std::cout << "ArcSim export has " << n_meshes << " meshes and " << n_obs << " obstacles." << std::endl;

	// If we're saving the frames out to a file (to make a video), we'll need a render window
	RenderWindow renderWindow;
	if( save_frames ){
		std::shared_ptr<DumpController> c = std::make_shared<DumpController>();
		renderWindow.set_controller( c );
		GLFWwindow* window = renderWindow.init();
		if( !window ){ return EXIT_FAILURE; }
	}

	// Load the obstacle meshes ahead of time. They are only transformed each frame.
	for( int i=0; i<n_obs; ++i ){
		obstacles.emplace_back( TriangleMesh::create() );

		std::string objfile = dir + "obs_" + pad_leading_zeros(2,i) + ".obj";
		if( !meshio::load_obj( obstacles.back().get(), objfile, false, false, false ) ){
			std::cerr << "\n**arcsimeToAlembic Error: Failed to load " << objfile << "\n" << std::endl;
		}
		std::string mesh_name = "obs"+std::to_string(i);

		obs_handles.emplace_back( exporter.add_object(mesh_name) );
		obs_aabb.emplace_back( obstacles.back()->bounds() );
	
		// Add the mesh to the render window
		if( save_frames ){
//			render_meshes.emplace_back( mcl::RenderMesh::create( obstacles.back() ) );
//			renderWindow.add_mesh( render_meshes.back(), RenderMesh::DYNAMIC );
		}
	}

	// Create the objects so that the alembic exporter knows about them
	for( int i=0; i<n_meshes; ++i ){
		std::string mesh_name = "mesh"+std::to_string(i);
		deform_handles.emplace_back( exporter.add_object(mesh_name) );
		meshes.emplace_back( TriangleMesh::create() );
		std::string meshstr = pad_leading_zeros(2,i);
		std::string objfile = dir + "000000_" + meshstr + ".obj";
		if( !meshio::load_obj( meshes.back().get(), objfile, false, false, false ) ){
			std::cerr << "\n**arcsimeToAlembic Error: Failed to load " << objfile << "\n" << std::endl;
			return EXIT_FAILURE;
		}

		// Add the mesh to the render window
		if( save_frames ){
			render_meshes.emplace_back( mcl::RenderMesh::create( meshes.back(), RenderMesh::DYNAMIC ) );
			renderWindow.add_mesh( render_meshes.back() );
		}
	}

	if( save_frames ){
		// Set the camera to a nice location based
		// on the meshes added to the renderWindow.
		renderWindow.nice_camera_location();

		// Game loop
		int frame_num = 0;
		bool stop = false;
		int n_render_meshes = render_meshes.size();
		int frame_output_count = 0;
		while( renderWindow.is_open() && !stop ){

			if( framedump_running ){
				std::cout << '\r';
				std::cout << "frame " << frame_num << std::flush;
				if( !load_frame( frame_num, stop ) ){
					return EXIT_FAILURE;
				}
				// Update rendering
				for( int i=0; i<n_render_meshes; ++i ){
					render_meshes[i]->load_buffers();
				}

				frame_num += skipframe;
				frame_output_count++;
			}

			renderWindow.draw();

			if( framedump_running ){
				// Save screenshot
				std::string ss_filename = make_screenshot_fn( frames_dir, frame_output_count );
				renderWindow.save_screenshot( ss_filename );
			}

			glfwPollEvents();

		} // end game loop


	} else {

		// Loop frames of the simulation and add them to the exporter
		std::cout << "Saving ArcSim frames to " << abc_name << std::endl;
		std::cout << "(This will take a while...)" << std::endl;
		bool stop = false;
		for( int frame_num=0; frame_num<999999; frame_num += skipframe ){
			std::cout << '\r';
			std::cout << "frame " << frame_num << std::flush;
			if( !load_frame( frame_num, stop ) ){
				return EXIT_FAILURE;
			}
			if( stop ){ break; } // all done!
		}
	}

	std::cout << "\nAll done!" << std::endl;

	return EXIT_SUCCESS;
}

void help(){

	std::cout << "\n====================\nUsage:" <<
		"\n\t -i <input directory>" <<
		"\n\t -dt <time step>" <<
		"\n\t -o <output file, optional>" <<
		"\n\t -dump <directory to save out frames, optional>" <<
		"\n\t -skip <frames to skip, optional>" <<
	"\n" << std::endl;
		

}

bool file_exists( const std::string &filename ){
	std::ifstream infile( filename.c_str() );
	if( !infile ){ return false; }
	return true;
}

int num_meshes( const std::string &dir ){
	int n_meshes = 0;
	for( int i=0; i<99; ++i ){
		std::stringstream obj;
		obj << dir << "000000_";
		if( i < 10 ){ obj << "0"; }
		obj << i << ".obj";
		if( file_exists( obj.str() ) ){
			n_meshes++;
		} else { break; }
	}
	if( n_meshes == 0 ){
		std::cerr << "\n**arcsimeToAlembic Error: Problem counting meshes\n" << std::endl;
	}
	return n_meshes;
}

int num_obs( const std::string &dir ){
	int n_obs = 0;
	for( int i=0; i<99; ++i ){
		std::stringstream obj;
		obj << dir << "obs_";
		if( i < 10 ){ obj << "0"; }
		obj << i << ".obj";
		if( file_exists( obj.str() ) ){
			n_obs++;
		} else { break; }
	}
	return n_obs;
}

std::string pad_leading_zeros( int num_digits, int num ){
	return std::string(num_digits-std::to_string(num).length(),'0') + std::to_string(num);
}

std::string get_name( const std::string &dir ){
	std::string filename = dir;
	filename.pop_back();
	std::size_t endir = filename.find_last_of("/\\");
	return filename.substr(endir+1);
}

bool test_conf( const std::string &dir ){
	std::stringstream conf_file;
	conf_file << dir << "conf.json";
	if( !file_exists( conf_file.str())  ){
		std::cerr << "\n**arcsimeToAlembic Error: No conf.json found, not an arcsim export\n" << std::endl;
		return false;
	}
	return true;
}
