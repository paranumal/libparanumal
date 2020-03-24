## libParanumal
An experimental set of finite element flow solvers for heterogeneous (GPU/CPU) systems. The initial development of libParanumal was performed by the [Parallel Numerical Algorithms Group at Virginia Tech](http://paranumal.com).   

libParanumal is funded in part by the US Department of Energy as part of the activities of the [Center for Efficient Exscale Discretizations](http://ceed.exascaleproject.org). 

Why is it called libParanumal ?: the high-order finite-element implementations in libParanumal are __spectrally__ accurate and rely heavily on __ghost__ elements for MPI communications.

If you use libParanumal as part of a research project see Section 8 below for papers to reference.

---
### 1. Overview 

Brief summary of major features:

A. Supported elements:
  - Triangles, quadrilaterals, tetrahedra, hexahedra.
  - Lagrange basis functions.
  
B. Mesh wrangling:
  - Gmsh format file loaders.
  - Load balanced geometric partitioning using space filling curves (Hilbert or Morton ordering). 
  
C. Time integrators:
  - Adaptive rate Dormand-Prince order 5 Runge-Kutta.
  - Low storage order 4 Runge-Kutta.
  - Single and Multirate Adams-Bashforth order 3.
  - Extrapolated Backwards Differencing order 3.
  
D. Iterative linear solvers:
  - Preconditioned (flexible) Conjugate Gradient method.
  - Non-blocking Preconditioned (flexible) Conjugate Gradient method.
  
E. Elliptic solver:
  - Linear Poisson and screened Poisson potential solvers.
  - GPU optimized matrix-vector products.
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
### 2. Code block diagram 
<img src="http://intranet.math.vt.edu/people/tcew/libPdiagramCrop.jpg" width="1024" >

---
### 3. OCCA dependency
`git clone https://github.com/libocca/occa`

#### 3-1. Build OCCA 
`cd occa`    
`export OCCA_DIR=${PWD}`  
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib`    
```make -j `nproc` ```    
`cd ../  `  

---
### 4. Clone: libParanumal
`git clone https://github.com/paranumal/libparanumal`

#### 4-1. Build all libParanumal solvers
`cd libparanumal`    
```make -j `nproc` ```

---
### 5. Running the codes: 

Each solver resides in its respective sub-directory in `solvers/`. Each solver sub-directory includes makefile, src directory, data directory (including header files for defining boundary conditions), okl kernel directory, and setups directory. The setups directory includes a number of example input files that specify input parameters for the solver.

#### 5-1. Build libParanumal elliptic solver
  
`cd libparanumal/solvers/elliptic`    
`make -j  `  

#### 5-2. Run elliptic example with provided quadrilateral set up file on a single device:
  
`./ellipticMain setups/setupQuad2D.rc`  

#### 5-3. Run the same example with two devices:

`mpiexec -n 2 ./ellipticMain setups/setupQuad2D.rc`  
 
---

### 6. License

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

### 7. References

Discontinuous Galerkin Boltzmann (bns) solver [arXiv version](https://arxiv.org/abs/1805.02082): `Karakus, A., Chalmers, N., Hesthaven, J.S. and Warburton, T., 2018. Discontinuous Galerkin Discretizations of the Boltzmann Equations in 2D: semi-analytic time stepping and absorbing boundary layers. arXiv preprint arXiv:1805.02082.`

Incompressible Navier-Stokes (discontinuous) Galerkin (ins) solver [publisher](https://doi.org/10.1016/j.jcp.2019.04.010), [arXiv version](https://arxiv.org/abs/1801.00246): `Karakus, A., Chalmers, N., Swirydowicz, K. and Warburton, T., 2017.A GPU accelerated discontinuous Galerkin incompressible flow solver, Journal of Computational Physics, Volume 390, Pages 380–404, 2019.`

Optimization of elliptic mat-vec operations for (elliptic) solver on hexes [publisher](https://doi.org/10.1177/1094342018816368), [arXiv version](https://arxiv.org/abs/1711.00903): `Kasia Swirydowicz, Noel Chalmers, Ali Karakus, and T. Warburton, Acceleration of tensor-product operations for high-order finite element methods, The International Journal of High Performance Computing Applications, Vol 33, Issue 4, 2019.`


