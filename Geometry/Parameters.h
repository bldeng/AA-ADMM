//  BSD 3-Clause License
//
//  Copyright (c) 2019, Bailin Deng
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
//  * Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Types.h"

class Parameters {
 public:

  Parameters()
      : iter(1),
        anderson_m(5),
        elasticity(sqrt(5000000.0)),
        time_step(0.033){
        //penalty_parameter(1){
  }

  virtual ~Parameters() {
  }

  // Parameters
  int iter;			// Number of iterations
  int anderson_m;     // Number of window size
  Scalar elasticity;  // Elasticity of spring
  Scalar time_step;
  //Scalar penalty_parameter;   // Penalty parameter for augmented Lagragian function, only used for geometry optimization


  // Load options from file
  bool load(const char* filename) {
    std::ifstream ifile(filename);
    if (!ifile.is_open()) {
      std::cerr << "Error while opening file " << filename << std::endl;
      return false;
    }

    std::string line;
    while (std::getline(ifile, line)) {
      std::string::size_type pos = line.find_first_not_of(' ');
      if (pos == std::string::npos) {
        continue;
      }

      // Check for comment line
      else if (line.at(pos) == '#') {
        continue;
      }

      std::string::size_type end_pos = line.find_first_of(' ');
      std::string option_str = line.substr(pos, end_pos - pos);
      std::string value_str = line.substr(end_pos + 1, std::string::npos);
      OptionInterpreter opt(option_str, value_str);

      load_option(opt);
    }

    std::cout << "Successfully loaded options from file " << filename
        << std::endl;

    return true;
  }

  void output() {
    std::cout << std::endl;
    std::cout << "====== Filter parameters =========" << std::endl;
    output_options();
    std::cout << "==================================" << std::endl;
    std::cout << std::endl;
  }

  // Check whether the parameter values are valid
  virtual bool valid_parameters() const {
    if (iter < 1) {
      std::cerr << "Error: Iterations must be at least 1" << std::endl;
      return false;
    }

    if (anderson_m < 0) {
      std::cerr << "Error: m must be non-negative" << std::endl;
      return false;
    }

    if (elasticity <= Scalar(0)) {
      std::cerr << "Error: elasticity must be positive" << std::endl;
      return false;
    }

    if (time_step <= Scalar(0)) {
      std::cerr << "Error: time step must be positive" << std::endl;
      return false;
    }

//    if (penalty_parameter <= Scalar(0)) {
//      std::cerr << "Error: penalty parameter must be positive" << std::endl;
//      return false;
//    }

    return true;
  }

 protected:

  class OptionInterpreter {
   public:
    OptionInterpreter(const std::string &option_str,
                      const std::string &value_str)
        : option_str_(option_str),
          value_str_(value_str) {
    }

    template<typename T>
    bool load(const std::string &target_option_name,
              T &target_option_value) const {
      if (option_str_ == target_option_name) {
        if (!load_value(value_str_, target_option_value)) {
          std::cerr << "Error loading option: " << target_option_name
              << std::endl;
          return false;
        }

        return true;
      } else {
        return false;
      }
    }

    template<typename EnumT>
    bool load_enum(const std::string &target_option_name, int enum_value_count,
                   EnumT &value) const {
      if (option_str_ == target_option_name) {
        int enum_int = 0;
        if (load_value(value_str_, enum_int)) {
          if (enum_int >= 0 && enum_int < enum_value_count) {
            value = static_cast<EnumT>(enum_int);
            return true;
          }
        }

        std::cerr << "Error loading option: " << target_option_name
            << std::endl;
        return false;
      } else {
        return false;
      }
    }

   private:
    std::string option_str_, value_str_;

    bool load_value(const std::string &str, double &value) const {
      try {
        value = std::stod(str);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Invalid argument: " << ia.what() << std::endl;
        return false;
      } catch (const std::out_of_range &oor) {
        std::cerr << "Out of Range error: " << oor.what() << std::endl;
        return false;
      }

      return true;
    }

    bool load_value(const std::string &str, float &value) const {
      try {
        value = std::stof(str);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Invalid argument: " << ia.what() << std::endl;
        return false;
      } catch (const std::out_of_range &oor) {
        std::cerr << "Out of Range error: " << oor.what() << std::endl;
        return false;
      }

      return true;
    }

    bool load_value(const std::string &str, int &value) const {
      try {
        value = std::stoi(str);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Invalid argument: " << ia.what() << std::endl;
        return false;
      } catch (const std::out_of_range &oor) {
        std::cerr << "Out of Range error: " << oor.what() << std::endl;
        return false;
      }

      return true;
    }

    bool load_value(const std::string &str, bool &value) const {
      int bool_value = 0;
      if (load_value(str, bool_value)) {
        value = (bool_value != 0);
        return true;
      } else {
        return false;
      }
    }
  };

  virtual bool load_option(const OptionInterpreter &opt) {
    return opt.load("Iterations", iter)
        || opt.load("AndersonM", anderson_m)
        || opt.load("SquareElasticity", elasticity)
        || opt.load("TimeStep", time_step);
        //|| opt.load("Penalty", penalty_parameter);
  }

  virtual void output_options() {
    std::cout << "Iteration number: " << iter << std::endl;

    std::cout << "Anderson m: " << anderson_m << std::endl;

    std::cout << "Square Elasticity (Only used in physics simulation): "
        << elasticity << std::endl;

    std::cout << "Time Step (Only used in physics simulation): " << time_step
        << std::endl;

    //std::cout << "Penalty parameter (only used in geometry optimiztion): " << penalty_parameter << std::endl;
  }

};

#endif // PARAMETERS_H
