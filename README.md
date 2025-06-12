# Rectangular Surface Parametrization

Code associated with the following publication:

[_Rectangular Surface Parametrization_](https://www.cs.cmu.edu/~kmcrane/Projects/RectangularSurfaceParameterization/RectangularSurfaceParameterization.pdf)\
Etienne Corman and Keenan Crane\
_Transaction on Graphics_, 2025

## Usage
The algorithm is launched with the script `run_RSP.m`. It will load an `.obj` file from the folder `Mesh/`. The output parametrization is exported to an `.obj` file in the folder `Results/`. 

## General options
The main script offers several options on the parametrization computation:
- `frame_field_type`: specifies how the initial cross field is computed. The parameter can take three values:
	- `'trivial'`: computes a cross field with fixed singularity indicies stored in `sing` using the algorithm described in *Trivial Connections on Discrete Surfaces*. Providing feasible singularities is the user's responsibility!
	- `'curvature'`: computes a curvature aligned cross field.
	- `'smooth'`: computes the smoothest frame field.
- `ifhardedge`: if `true` hard-edges are detected and used as field/parametrization alignment constraints.
- `ifboundary`: if `true` the field/parametrization is constrained to align with the surface boundary. Boundary alignment is necessary for the quantization step.
- `ifseamless_const`: if `true` the parametrization is computed with strong seamlessness constraints. This is necessary for the quantization step.
- `ifquantization`: if `true` the parametrization is quantized in post-processing. The quantization requires boundary alignment (ie `ifboundary = true`) and exact seamless map as input (ie `ifseamless_const = true`).

## Objective functions
Various objective functions (see Section 5.1.1) can be chosen using the variable `energy_type` whose options are accessible with the structure `weight`. The variable `energy_type` takes three values:
- `'distortion'`: quadratic energy on `u` and `v` controlled by the real value `weight.w_conf_ar`:
	- `weight.w_conf_ar = 0`: as *area-preserving* as possible;
	- `weight.w_conf_ar = 0.5`: as *isometric* as possible;
	- `weight.w_conf_ar = 1`: as *conformal* as possible.
- `'chebyshev'`: promotes Chebyshev net (see Section 6.6 of the paper).
- `'alignment'`: penalizes the frame angle `ang` to be close to the target angle `weight.ang_dir`. The aspect ratio `v` can also be closed to the target `weight.aspect_ratio`. The weight of each energy is set by `weight.w_ang` and `weight.w_ratio` respectively.

The variable `weight.w_gradv` sets the weight of the regularizer (see Section 5.1.2).

## Examples
Three use case are provided:
1. Smoothest cross field with hard-edge constraints:
```
mesh_name = 'B36';
frame_field_type = 'smooth';
ifhardedge   = true;
ifboundary   = true;
ifseamless_const = true;
ifquantization = true;
energy_type = 'distortion';
``` 
2. Curvature aligned parametrization:
```
mesh_name = 'pig';
frame_field_type = 'curvature'
ifhardedge   = false;
ifboundary   = true;
ifseamless_const = true;
ifquantization = true;
energy_type = 'alignment';
``` 
3. Chebyshev net computation with boundary alignment:
```
mesh_name = 'SquareMyles';
frame_field_type = 'trivial'
ifhardedge   = false;
ifboundary   = true;
ifseamless_const = true;
ifquantization = false;
energy_type = 'chebyshev';
``` 

## Quantization
The quantization step turn a seamless map into an *integer* seamless map. This is done using *Quad Mesh Quantization Without a T-Mesh* by Yoann Coudert-Osmont, David Desobry, Martin Heistermann, David Bommes, Nicolas Ray and Dmitry Sokolov. 

The source code is located in the folder `QuantizationYoann/` and is gracefully provided by [Yoann Coudert-Osmont](https://perso.ens-lyon.fr/yoann.coudert-osmont/). It must be compiled before use.

#### Compilation
The compilation requires CMake and a C++ compiler. From the folder `QuantizationYoann/`, run in a terminal:
```
mkdir build/
cd build/
cmake ..
make
```
The program `Quantization` is generated in the `build/` subdirectory.

#### Limitations
The quantization step will fail when one of these three situations occur:
1. the parametrization is not aligned with boundary;
2. the parametrization is not *exactly* seamless;
3. the parametrization is not locally injective (ie some Jacobian matrix have negative determinant).

Keep in mind that the algorithm does not provide any guaranties regarding local injectivity of rectangular parametrizations.

## Quad-mesh extraction
A quad-mesh can be extracted directly from an integer seamless map. We do not provide code for this step. The meshes shown in the paper were obtained with [libQEx – A Robust Quad Mesh Extractor](https://github.com/hcebke/libQEx).

## License
Copyright (C) 2025, Etienne CORMAN and Keenan CRANE \
SPDX-License-Identifier: AGPL-3.0-or-later

If you make use of this code in scientific work we kindly ask you to cite our paper:
```
@article{Corman:2025:RSP,
author = {Corman, Etienne and Crane, Keenan},
title = {Rectangular Surface Parameterization},
journal = {ACM Trans. Graph.},
volume = {44},
number = {4},
year = {2025},
publisher = {ACM},
address = {New York, NY, USA},
}
``` 

*Commercial licensing* under negotiable terms is available upon request. Please send an email to etienne.corman@cnrs.fr and kmcrane@cs.cmu.edu if you are interested.