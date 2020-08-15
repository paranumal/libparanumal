![Build](https://github.com/paranumal/libparanumal/workflows/Build/badge.svg)
[![codecov](https://codecov.io/gh/paranumal/libparanumal/branch/master/graph/badge.svg)](https://codecov.io/gh/paranumal/libparanumal)

## libParanumal
An experimental set of finite element flow solvers for heterogeneous (GPU/CPU) systems. The initial development of libParanumal was performed by the [Parallel Numerical Algorithms Group at Virginia Tech](http://paranumal.com).

libParanumal is funded in part by the US Department of Energy as part of the activities of the [Center for Efficient Exscale Discretizations](http://ceed.exascaleproject.org).

Why is it called libParanumal?: the high-order finite-element implementations in libParanumal are __spectrally__ accurate and rely heavily on __ghost__ elements for MPI communications.

### 1. How to cite the libParanumal project:
If you use any part of libParanumal in your research project including variants on the included compute kernels please use the following bibliographical reference:

<pre>
@MISC{ChalmersKarakusAustinSwirydowiczWarburton2020,
      author = "Chalmers, N. and Karakus, A. and Austin, A. P. and Swirydowicz, K. and Warburton, T.",
      title = "{libParanumal}: a performance portable high-order finite element library",
      year = "2020",
      url = "https://github.com/paranumal/libparanumal",
      note = "Release 0.3.1"
      }
</pre>

see the [references](#10-references) section below for additional papers to reference about various aspects of the library.

---
### 2. How you can help out:
libParanumal is a community project. Please help improve the library by submitted an [Issue](https://github.com/paranumal/libparanumal/issues)  if you notice any unexpected behavior, discover a bug, have problems installing/running/trouble shooting your installation. It benefits us as a community when issues and feature requests are shared so that we can understand how to improve the library.

Please submit feature requests as an [Issue](https://github.com/paranumal/libparanumal/issues) for the consideration of all contributors. Likewise if you wish to submit a code change please make a GitHub [Pull Request](https://github.com/paranumal/libparanumal/pulls).

---
### 3. Overview

Brief summary of major features:

A. Supported elements:
  - Triangles, quadrilaterals, tetrahedra, hexahedra.
  - Lagrange basis functions.

B. Mesh wrangling:
  - Gmsh format file loaders.
  - Load balanced geometric partitioning using space filling curves (Hilbert or Morton ordering).

C. Time integrators:
  - Adaptive rate Dormand-Prince order 5 Runge-Kutta.
  - Low storage explicit Runge-Kutta order 4.
  - Single and Multirate Adams-Bashforth order 3.
  - Extrapolated Backwards Differencing order 3.

D. Iterative linear solvers:
  - Preconditioned (flexible) Conjugate Gradient method.
  - Non-blocking Preconditioned (flexible) Conjugate Gradient method.

E. Elliptic solver:
  - Linear Poisson and screened Poisson potential solvers.
  - GPU-optimized matrix-vector products.
  - p-type multigrid, algebraic multigrid, low-order SEMFEM, and Jacobi preconditioning.
  - Matrix-free p-multigrid for fine levels of multigrid hierarchy.

F. Heterogeneous accelerated flow solvers:
  - Compressible Navier-Stokes solver with:
     * Upwind discontinuous Galerkin discretization in space.
     * Optional isothermal equation of state.
  - Galerkin-Boltzmann gas dynamics solver with:
     * Penalty flux DG discretization in space.
     * Adaptive semi-analytic (pointwise exponential) integration in time.
     * Multiaxial quasi-perfectly matched absorbing layer far field boundary condition.
  - Incompressible Navier-Stokes solver with:
     * Choice of continuous FEM or interior penalty DG in space.
     * Extrapolation-BDF integration in time.
     * Sub-cycling (Operator Integration Factor Splitting) for advection.

G. Dependencies:
   - Message Passing Interface (MPI v3.0 or higher).
      * The libParanumal makefiles assume that mpic++ is installed and visible in your path.
   - Open Concurrent Compute Abstraction (OCCA)
      * OCCA must be installed.
      * OCCA will try to detect if any of these execution models are installed: OpenMP, CUDA, OpenCL, and/or HIP.
      * By default, if OCCA does not detect a chosen mode of execution it will default to Serial execution.
      * You will need to adjust the libParnumal setup input files to choose the execution model and compute device appropriate for your system.
      * The OCCA github repo is [here](https://github.com/libocca/occa)
      * The OCCA webpage is [here](http://libocca.org)


---
### 4. Code block diagram
<img src="http://intranet.math.vt.edu/people/tcew/libPdiagramCrop.jpg" width="1024" >

---
### 5. OCCA dependency
`git clone https://github.com/libocca/occa`

#### 5-1. Build OCCA
`cd occa`
`export OCCA_DIR=${PWD}`
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib`
```make -j `nproc` ```
`cd ../  `

---
### 6. Required Libraries
libParanumal requires installed BLAS and LAPACK libraries. By default, the build system will look for `libblas` and `liblapack` in your default library search paths. The library paths can also be manually specified in `make.top` with the `LIBP_BLAS_DIR` and `LIBP_LAPACK_DIR` variables.

Some Linux distributions will package BLAS and LAPACK libraries. For example, on Ubuntu systems these libraries can be installed via
```sudo apt install libblas-dev liblapack-dev```

---
### 7. Clone: libParanumal
`git clone https://github.com/paranumal/libparanumal`

#### 7-1. Build all libParanumal solvers
`cd libparanumal`
```make -j `nproc` ```

---
### 8. Running the codes:

Each solver resides in its respective sub-directory in `solvers/`. Each solver sub-directory includes makefile, src directory, data directory (including header files for defining boundary conditions), okl kernel directory, and setups directory. The setups directory includes a number of example input files that specify input parameters for the solver.

#### 8-1. Build libParanumal elliptic solver

`cd libparanumal/solvers/elliptic`
`make -j  `

#### 8-2. Run elliptic example with provided quadrilateral set up file on a single device:

`./ellipticMain setups/setupQuad2D.rc`

#### 8-3. Run the same example with four devices:

`mpiexec -n 4 ./ellipticMain setups/setupQuad2D.rc`

---

### 9. License

The MIT License (MIT)

Copyright (c) 2017-2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### 10. References

Discontinuous Galerkin Boltzmann (bns) solver: [publisher](https://doi.org/10.1016/j.jcp.2019.03.050), [arXiv version](https://arxiv.org/abs/1805.02082): `Karakus, A., Chalmers, N., Hesthaven, J.S. and Warburton, T., 2018. Discontinuous Galerkin Discretizations of the Boltzmann Equations in 2D: semi-analytic time stepping and absorbing boundary layers, Journal of Computational Physics, Volume 390, Pages 175–202.`

Incompressible Navier-Stokes (discontinuous) Galerkin (ins) solver: [publisher](https://doi.org/10.1016/j.jcp.2019.04.010), [arXiv version](https://arxiv.org/abs/1801.00246): `Karakus, A., Chalmers, N., Swirydowicz, K. and Warburton, T., 2017.A GPU accelerated discontinuous Galerkin incompressible flow solver, Journal of Computational Physics, Volume 390, Pages 380–404, 2019.`

Optimization of elliptic mat-vec operations for (elliptic) solver on hexes: [publisher](https://doi.org/10.1177/1094342018816368), [arXiv version](https://arxiv.org/abs/1711.00903): `Swirydowicz, K., Chalmers, N., Karakus, A., and Warburton, T. 2019. Acceleration of tensor-product operations for high-order finite element methods, The International Journal of High Performance Computing Applications, Vol 33, Issue 4.`

Low-order preconditioning of triangular elements (elliptic precon): [publisher](https://epubs.siam.org/doi/abs/10.1137/17M1149444): `Chalmers, N. and Warburton, T. 2018. Low-order preconditioning of high-order triangular finite elements, SIAM Journal on Scientific Computing, Vol 40, Issue 6, Pages A4040-A4059`
