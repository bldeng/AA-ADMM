## AA-ADMM: Anderson Acceleration for ADMM


This is the source code for the algorithms in the following paper:

* Juyong Zhang<sup>\*</sup>, Yue Peng<sup>\*</sup>, Wenqing Ouyang<sup>\*</sup>, Bailin Deng. 2019. Accelerating ADMM for Efficient Simulation and Optimization. arXiv:1909.00470. (<sup>\*</sup>joint first authors).

### Content

The folders contain the following code:

* `admm_anderson_hard_zxu` and `admm_ander_xzu`: accelerated solver for physics simulation, based on the [implementation of \[Overby et al. 2017\]](https://github.com/mattoverby/admm-elastic).

* Geometry: accelerated ADMM for wire mesh optimization and planar quad mesh optimization, based on the formulation of [\[Deng et al. 2015\]](http://orca.cf.ac.uk/98569/1/98569Interactive%20design%20exploration%20for%20constrained%20meshes.pdf).

Please checkout the README file in each folder for instructions on usage.


### License

The code is released under BSD 3-Clause License.

### Contact

For comments and questions, please contact Yue Peng <<echoyue@mail.ustc.edu.cn>> or Bailin Deng <<bldeng@gmail.com>>.
