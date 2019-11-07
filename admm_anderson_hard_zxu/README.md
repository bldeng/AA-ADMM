##Anderson Acceleration for ADMM

### Compiling

The Code has been tested on the following. Systems and compilers:

* Ubuntu 18.04 using gcc 7.3. And it needs MKL, or it may slow down the performance.

Follow the following steps to compile the code:

* Create a folder `build` within the root directories of the code

* Run cmake to generate the build files inside the build folder, and compile the source code:
	* On linux or mac, run the following commands within the build folder:

    	```
    	$ cmake -DCMAKE_BUILD_TYPE=Release ..
    	$ make
    	```
	* On mac, the default Apple clang compiler does not support OpenMP. To enable OpenMP, first install the homebrew GCC:

		```
		$ brew install gcc
		```
	  This should install gcc 8, with compiler commands gcc-8 and g++8. Then run the following command within the build folder:

	  	```
	  	$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 ..
	  	$ make
	  	```
	    

* Afterwards, there should be a folder `Asia2019` generated under `admm_anderson/build/samples/`, and there are few executable files in it.

	* Create new result folder, output message will be generated here:

		```
       $ mkdir ./samples/Asia2019/result
       ```

### Usage
* Some examples:

	```
	$ ./samples/Asia2019/beams -a 0
	$ ./samples/Asia2019/beams -am 3
	$ ./samples/Asia2019/plinkohit -am 6
	```

	* -a 0: not accelerate.
	* -am 1: use anderson acceleration and set window size m=1

* Simulation demo instructions:
	* Press 'p' to simulate one frame.
	* Press 'space' to start demo.

* A script for test has been included. Here are some instructions:
	* Put both `testAndersonADMM` and `testParam.txt` into the folder `Asia2019`.
	* Run 
	
		```
		$ cd ./samples/Asia2019/
		$ ./testAndersonADMM
		```
   * Press 'p' and close the window after it finishs simiulation.
   * Repeat the last step.
   * You can find results in the subfolder `result`.
