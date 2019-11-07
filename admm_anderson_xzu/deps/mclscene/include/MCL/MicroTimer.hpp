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

#ifndef MCL_MICROTIMER_H
#define MCL_MICROTIMER_H 1

#include <chrono>
#include <unordered_map>

// Example:
//
//	mcl::MicroTimer timer;
//
//	... do work ...
//
//	double elapsed_microseconds = timer.elapsed_us();
//	double elapsed_milliseconds = timer.elapsed_ms();
//	double elapsed_seconds = timer.elapsed_s();
//
//	Notes:
//		You may want to change the clock type based on application
//
namespace mcl {

// I call it MicroTimer to avoid name clashes.
// Also because it can do microseconds, though you'd need the hrc
// instead of stead.
class MicroTimer {
	//typedef std::chrono::high_resolution_clock C;
	typedef std::chrono::steady_clock C;
	typedef double T;
public:

	MicroTimer() : start_time( C::now() ){}

	// Resets the timer
	void reset() { start_time = C::now(); }

	// Return time elapsed in seconds
	T elapsed_s() const {
		curr_time = C::now();
		std::chrono::duration<T> durr = curr_time-start_time;
		return durr.count();
	}

	// Return time elapsed in milliseconds
	T elapsed_ms() const {
		curr_time = C::now();
		std::chrono::duration<T, std::milli> durr = curr_time-start_time;
		return durr.count();
	}

	// Return time elapsed in microseconds
	T elapsed_us() const {
		curr_time = C::now();
		std::chrono::duration<T, std::micro> durr = curr_time-start_time;
		return durr.count();
	}

private:
	std::chrono::time_point<C> start_time;
	mutable std::chrono::time_point<C> curr_time;

}; // end class MicroTimer


//
// TimerManager is just a collection of named timers.
// Eventually I'll add more functionality as needed.
//
class TimerManager {
public:
	// Default unit is millesecond
	TimerManager() : unit(1) {}

	// All timers are the same unit ("s", "ms", or "us");
	inline void set_unit( const std::string &unit_ ){
		if( unit_=="s"){ unit=0; }
		else if( unit_=="s" ){ unit=1; }
		else if( unit_=="us" ){ unit=2; }
	}

	// Begin a named timer
	inline void start( const std::string &name ){
		if( totals.count(name) == 0 ){
			totals[name] = 0.0;
			counts[name] = 0;
		}
		current[name] = MicroTimer();
	}

	// Stop a named timer
	inline void stop( const std::string &name ){
		std::unordered_map<std::string, MicroTimer>::iterator it = current.find(name);
		if( it == current.end() ){ return; } // if the timer was never started...
		double update = 0.0; // find current elapsed time
		switch(unit){
			case 0:{ update = it->second.elapsed_s(); } break;
			case 1:{ update = it->second.elapsed_ms(); } break;
			case 2:{ update = it->second.elapsed_us(); } break;
		}
		current.erase(it); // remove the ongoing timer
		totals[name] += update;
		counts[name] += 1;
	}

	// Adds to a named timer
	inline void add( const std::string &name, double t ){
		if( totals.count(name) == 0 ){
			totals[name] = 0.0;
			counts[name] = 0;
		}
		totals[name] += t;
		counts[name] += 1;
	}

	// Get total elapsed time for a named timer
	inline double total( const std::string &name ){
		std::unordered_map<std::string, double>::iterator it = totals.find(name);
		if( it == totals.end() ){ return 0.0; } // if the timer was never started...
		return it->second;
	}

	// Get average elapsed time for a named timer
	inline double average( const std::string &name ){
		std::unordered_map<std::string, double>::iterator it = totals.find(name);
		if( it == totals.end() ){ return 0.0; } // if the timer was never started...
		double count = (double)counts[it->first];
		return it->second / count;
	}

	// Prints the average of all named timers
	inline void print_averages(){
		std::unordered_map<std::string, double>::iterator it = totals.begin();
		for( ; it != totals.end(); ++it ){
			// Should never fail, totals and counts set at same time...
			double count = (double)counts[it->first];
			double avg = it->second / count;
			std::cout << "Avg time " << it->first << ": " << avg;
			switch(unit){
				case 0:{ std::cout << " s"; } break;
				case 1:{ std::cout << " ms"; } break;
				case 2:{ std::cout << " us"; } break;
			}
			std::cout << std::endl;
		} // end loop named times
	}

private:
	int unit; // 0=s, 1=ms, 2=us
	std::unordered_map<std::string, double> totals;	
	std::unordered_map<std::string, int> counts;
	std::unordered_map<std::string, MicroTimer> current;
};


} // end namespace mcl


#endif
