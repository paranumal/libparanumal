## List of examples

The libparanumal solver subdirectories include finite element discretizations of the following PDEs:

**acoustics**
* Linearized Euler equations:
  * 2D:
       * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}$$
       * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}$$
       * $$\frac{\partial p}{\partial t} = -\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}$$
  * 3D:
       * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}$$
       * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}$$
       * $$\frac{\partial w}{\partial t} = -\frac{\partial p}{\partial z}$$
       * $$\frac{\partial p}{\partial t} = -\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}-\frac{\partial w}{\partial z}$$

**elliptic**
* Screened Poisson potential problem
  * 2D:
    * PDE: $$\lambda u - \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)=f \;(x,y,z)\in\Omega$$
    * Dirichlet BC: $$u=g(x,y) \mbox{ for } (x,y)\in\partial\Omega^D$
    * Neumann BC:   $$n_x\frac{\partial u}{\partial x}+n_y\frac{\partial u}{\partial y}=h(x,y) \mbox{  for } (x,y)\in\partial\Omega^N$$

  * 3D:
    * PDE: $$\lambda u - \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right)=f$$
    * Dirichlet BC: $$u=g(x,y,z) \mbox{ for } (x,y,z) \in \partial\Omega^D$$
    * Neumann BC:   $$n_x\frac{\partial u}{\partial x}+n_y\frac{\partial u}{\partial y}+n_z\frac{\partial u}{\partial z}=h(x,y) \mbox{ for } (x,y) \in \partial\Omega^N$$
