# SampleDataCBF: Sampled-Data Control Barrier Function
This repository contains source code for the paper [Control barrier function meets interval analysis: Safety-critical control with measurement and actuation uncertainties](https://ieeexplore.ieee.org/document/9867681). The arxiv version of the paper can be found [here](http://128.84.4.34/abs/2110.00915).
The repository contains MATLAB code that has been tested in MATLAB R2022b.

If you use this code, we would appreciate it if you cited the paper as follows:
```
@inproceedings{zhang2022control,
  title={Control barrier function meets interval analysis: Safety-critical control with measurement and actuation uncertainties},
  author={Zhang, Yuhao and Walters, Sequoyah and Xu, Xiangru},
  booktitle={2022 American Control Conference (ACC)},
  pages={3814--3819},
  year={2022},
  organization={IEEE}
}
```

For more information about our work, please visit [ARC Lab@UW-Madison](https://xu.me.wisc.edu/).

## MATLAB
### Dependencies & Installation
The MATLAB implementation requires a working installation of
[CORA](https://tumcps.github.io/CORA/). This repository also contains a modified fork of [SparsePOP](https://sparsepop.sourceforge.io/) in the folder `SparsePOP303`. We recommend using [MOSEK](https://www.mosek.com/) as the Semi-Definite Programming (SDP) solver for [SparsePOP](https://sparsepop.sourceforge.io/), noting that you will need a license (free academic license available). However, other SDP solvers supported by [SparsePOP](https://sparsepop.sourceforge.io/) should also work, such as [SeDuMi](https://github.com/sqlp/sedumi) and [SDPT3](https://github.com/sqlp/sdpt3).

You might also find the following MATHWORKS toolboxes useful/necessary:
* Optimization Toolbox

The package is lightweight and there is no installation beyond adding the following folders to
your path:
```
addpath(genpath(userpath+"SampleDataCBF"))
```
Note that if you want to use [MOSEK](https://www.mosek.com/), you will also need add the following folders to
your path:
```
addpath(genpath(userpath+"Mosek/10.0/toolbox"))
```

### Example

For a quick-start, inspect the file `sampleCBF_sim_main.m` and run it.
You should see the following plot.

<p align="center">
  <img src="https://github.com/wisc-arclab/SampleDataCBF/blob/master/img/example_jankovic.png" width="500" alt="Jancovic Example Path">
</p>

You can also test the mass-spring-damper example by setting
```
param.example = 'spring';
```
And you will get the following plots.

<p float="center">
  <img src="https://github.com/wisc-arclab/SampleDataCBF/blob/master/img/example_spring_states.png" width="400" alt="Spring Example States">
  <img src="https://github.com/wisc-arclab/SampleDataCBF/blob/master/img/example_spring_x1.png" width="400" alt="Spring Example x1">
</p>

The simulation of the linearized 6-dimension quadcopter model can be executed by running the file `sampleCBF_sim_DI_main.m` and you will get the following plots.

<p float="center">
  <img src="https://github.com/wisc-arclab/SampleDataCBF/blob/master/img/example_di.png" width="400" alt="DI Example States">
  <img src="https://github.com/wisc-arclab/SampleDataCBF/blob/master/img/example_di_y.png" width="400" alt="DI Example y">
</p>


## Acknowledgements
This work was supported by the University of Wisconsin.