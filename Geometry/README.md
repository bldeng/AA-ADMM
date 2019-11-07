##Anderson Acceleration for ADMM


### Compiling

The Code has been tested on Ubuntu 18.04 using gcc 7.3. To compile the code:

* Create a folder “build” within the root directories of the code

	* Run cmake to generate the build files inside the build folder, and compile the source code:

   * On linux or mac, run the following commands within the build folder:
    	$ cmake -DCMAKE_BUILD_TYPE=Release ..
    	$ make
    	
	* On mac, the default Apple clang compiler does not support OpenMP. To enable OpenMP, first install the homebrew GCC:
	
		```
		$ brew install gcc
		```
	
	  This should install gcc 8, with compiler commands gcc-8 and g++8. Then run the following command within the build folder:
	  
	  ```
	    $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 ..  
	    $ make
	    ```
	    
* Afterwards, there are two exectuables in the build folder: `WireMeshOpt` and `PlanarityOpt`.

### Usage
Use the commands
	```
	$ WireMeshOpt INPUT_POLY_MESH  REF_TRI_MESH  OPTIONS_FILE  OUTPUT_MESH
	```
or
	```
	$ PlanarityOpt INPUT_POLY_MESH  REF_TRI_MESH  OPTIONS_FILE  OUTPUT_MESH
	```
where

* `INPUT_POLY_MESH ` is the input polygonal mesh.
*  `REF_TRI_MESH ` is reference triangle mesh.
*  `OPTIONS_FILE ` is the option file (see `Options.txt` for an example).
*  `OUTPUT_MESH ` is a file to save the optimized output mesh.

You can find some test models in the folder `Geometry_model`.
