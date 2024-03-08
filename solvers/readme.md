## List of examples

The libparanumal solver subdirectories include finite element discretizations of the following PDEs:

**acoustics**
* Time dependent linearized Euler equations:
  * 2D:
       * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}$$
       * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}$$
       * $$\frac{\partial p}{\partial t} = -\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}$$
  * 3D:
       * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}$$
       * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}$$
       * $$\frac{\partial w}{\partial t} = -\frac{\partial p}{\partial z}$$
       * $$\frac{\partial p}{\partial t} = -\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}-\frac{\partial w}{\partial z}$$

**advection**
* Linear advection equation:
  * 2D:
       * $$\frac{\partial u}{\partial t} = v_x\frac{\partial u}{\partial x}+v_y\frac{\partial u}{\partial y}$$
  * 3D:
       * $$\frac{\partial u}{\partial t} = v_x\frac{\partial u}{\partial x}+v_y\frac{\partial u}{\partial y}+v_z\frac{\partial u}{\partial z}$$

**bns**
* Boltzmann Navier-Stokes equations of gas dynamics
  * 2D:
  * PDE: 
    $$\frac{\partial q_0}{\partial t} = -c\left(\frac{\partial q_1}{\partial x} + \frac{\partial q_2}{\partial y}\right)$$
    $$\frac{\partial q_1}{\partial t} = -c\left(\frac{\partial q_0}{\partial x} + \sqrt{2}\frac{\partial q_4}{\partial x} + \frac{\partial q_3}{\partial y}\right)$$
    $$\frac{\partial q_2}{\partial t} = -c\left(\frac{\partial q_3}{\partial x} + \frac{\partial q_0}{\partial y} + \sqrt{2}\frac{\partial q_5}{\partial y}\right)$$
    $$\frac{\partial q_3}{\partial t} = -c\left(\frac{\partial q_2}{\partial x} + \frac{\partial q_1}{\partial y}\right)-\frac{1}{\tau}\left(\frac{q_3-q_1q_2}{q_0}\right)$$
    $$\frac{\partial q_4}{\partial t} = -c\sqrt{2}\left(\frac{\partial q_1}{\partial x}\right)-\frac{1}{\tau}\left(\frac{q_4-q_1^2}{q_0\sqrt{2}}\right)$$
    $$\frac{\partial q_5}{\partial t} = -c\sqrt{2}\left(\frac{\partial q_2}{\partial y}\right)-\frac{1}{\tau}\left(\frac{q_5-q_2^2}{q_0\sqrt{2}}\right)$$    

**cns**
* Compressible Navier-Stokes (does not include: limiting, artificial viscosity, entropy stable formulation)
  * 2D:
    $$\frac{\partial \rho}{\partial t} = -\frac{\partial\rho{}u}{\partial x}  -\frac{\partial\rho{}v}{\partial y}$$
    $$\frac{\partial \rho{}u}{\partial t} = -\frac{\partial}{\partial x}\left(\rho{}u^2+p-\mu\tau_{11}\right)  -\frac{\partial}{\partial y}\left(\rho{}uv-\mu\tau_{12}\right)$$
    $$\frac{\partial \rho{}v}{\partial t} = -\frac{\partial}{\partial x}\left(\rho{}uv-\mu\tau_{21}\right)  -\frac{\partial}{\partial y}\left(\rho{}v^2+p-\mu\tau_{22}\right)$$
     $$\frac{\partial E}{\partial t} =
    -\frac{\partial}{\partial x}\left(u(E+p)-\mu\left(u\tau_{11}+v\tau_{12}\right)\right)
    -\frac{\partial}{\partial y}\left(v(E+p)-\mu\left(u\tau_{12}+v\tau_{22}\right)\right)$$
    where:
    $$\tau_{11} = 2\frac{\partial u}{\partial x} - \frac{2}{3}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)$$
    $$\tau_{22} = 2\frac{\partial v}{\partial y} - \frac{2}{3}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)$$
    $$\tau_{12} = \tau_{21} = \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}$$
    $$E = \frac{p}{\gamma-1} + \frac{1}{2}\rho\left(u^2+v^2\right)$$
 * 3D:
    $$\frac{\partial \rho}{\partial t} = -\frac{\partial\rho{}u}{\partial x}  -\frac{\partial\rho{}v}{\partial y} -\frac{\partial\rho{}w}{\partial z}$$
    $$\frac{\partial \rho{}u}{\partial t} =
   -\frac{\partial}{\partial x}\left(\rho{}u^2+p-\mu\tau_{11}\right)
   -\frac{\partial}{\partial y}\left(\rho{}uv-\mu\tau_{12}\right)
    -\frac{\partial}{\partial z}\left(\rho{}uw-\mu\tau_{13}\right)$$
    $$\frac{\partial \rho{}v}{\partial t} =
   -\frac{\partial}{\partial x}\left(\rho{}uv-\mu\tau_{21}\right)
   -\frac{\partial}{\partial y}\left(\rho{}v^2+p-\mu\tau_{22}\right)
   -\frac{\partial}{\partial z}\left(\rho{}wv-\mu\tau_{23}\right)$$
       $$\frac{\partial \rho{}w}{\partial t} =
   -\frac{\partial}{\partial x}\left(\rho{}uw-\mu\tau_{31}\right)
   -\frac{\partial}{\partial y}\left(\rho{}vw-\mu\tau_{32}\right)
   -\frac{\partial}{\partial z}\left(\rho{}w^2+p-\mu\tau_{33}\right)$$
     $$\frac{\partial E}{\partial t} =
    -\frac{\partial}{\partial x}\left(u(E+p)-\mu\left(u\tau_{11}+v\tau_{12}+w\tau_{13}\right)\right)
    -\frac{\partial}{\partial y}\left(v(E+p)-\mu\left(u\tau_{21}+v\tau_{22}+w\tau_{23}\right)\right)
   -\frac{\partial}{\partial z}\left(w(E+p)-\mu\left(u\tau_{31}+v\tau_{32}+3\tau_{33}\right)\right)$$
    where:
    $$\tau_{11} = 2\frac{\partial u}{\partial x} - \frac{2}{3}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}\right)$$
    $$\tau_{22} = 2\frac{\partial v}{\partial y} - \frac{2}{3}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}\right)$$
    $$\tau_{33} = 2\frac{\partial w}{\partial z} - \frac{2}{3}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}\right)$$
    $$\tau_{12} = \tau_{21} = \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}$$
    $$\tau_{13} = \tau_{31} = \frac{\partial u}{\partial z}+\frac{\partial w}{\partial x}$$
    $$\tau_{23} = \tau_{32} = \frac{\partial v}{\partial z}+\frac{\partial w}{\partial y}$$
    $$E = \frac{p}{\gamma-1} + \frac{1}{2}\rho\left(u^2+v^2+w^2\right)$$
    
**elliptic**
* Screened Poisson potential problem
  * 2D:
    * PDE: $$\lambda u - \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)=f \mbox{ for } (x,y,z)\in\Omega$$
    * Dirichlet BC: $$u=g(x,y) \mbox{ for } (x,y)\in\partial\Omega^D$$
    * Neumann BC:   $$n_x\frac{\partial u}{\partial x}+n_y\frac{\partial u}{\partial y}=h(x,y) \mbox{  for } (x,y)\in\partial\Omega^N$$

  * 3D:
    * PDE: $$\lambda u - \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right)=f$$
    * Dirichlet BC: $$u=g(x,y,z) \mbox{ for } (x,y,z) \in \partial\Omega^D$$
    * Neumann BC:   $$n_x\frac{\partial u}{\partial x}+n_y\frac{\partial u}{\partial y}+n_z\frac{\partial u}{\partial z}=h(x,y) \mbox{ for } (x,y) \in \partial\Omega^N$$

**fokkerPlanck**
* Time-domain advection diffusion modeled with the Fokker-Planck equations:
  * 2D:
     * $$\frac{\partial u}{\partial t} = v_x \frac{\partial u}{\partial x} + v_y \frac{\partial u}{\partial y} + \mu \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)$$
  * 3D:
     * $$\frac{\partial u}{\partial t} = v_x \frac{\partial u}{\partial x} + v_y \frac{\partial u}{\partial y} + v_z \frac{\partial u}{\partial z} + \mu \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right)$$

**ins**
* Time-dependent incompressible Navier-Stokes equations:
  * 2D:
     * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}-u \frac{\partial u}{\partial x} - v \frac{\partial u}{\partial y} + \mu \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)$$
     * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}-u \frac{\partial v}{\partial x} - v \frac{\partial v}{\partial y} + \mu \left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)$$
     * $$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}=0$$
  * 3D:
     * $$\frac{\partial u}{\partial t} = -\frac{\partial p}{\partial x}-u \frac{\partial u}{\partial x} - v \frac{\partial u}{\partial y}- w \frac{\partial u}{\partial z} + \mu \left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right)$$
     * $$\frac{\partial v}{\partial t} = -\frac{\partial p}{\partial y}-u \frac{\partial v}{\partial x} - v \frac{\partial v}{\partial y}- w \frac{\partial u}{\partial z} + \mu \left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}+\frac{\partial^2 v}{\partial z^2}\right)$$
     *  $$\frac{\partial w}{\partial t} = -\frac{\partial p}{\partial z}-u \frac{\partial w}{\partial x} - v \frac{\partial w}{\partial y}- w \frac{\partial w}{\partial z} + \mu \left(\frac{\partial^2 w}{\partial x^2}+\frac{\partial^2 w}{\partial y^2}+\frac{\partial^2 w}{\partial z^2}\right)$$
     * $$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}+ \frac{\partial w}{\partial z}=0$$
   
**maxwell**
* Time-domain Maxwell equations of electromagnetics:
  * 2D (transverse mode):
       * $$\frac{\partial H_x}{\partial t} = -\frac{\partial E_z}{\partial y}$$
       * $$\frac{\partial H_y}{\partial t} =  \frac{\partial E_z}{\partial x}$$
       * $$\frac{\partial E_z}{\partial t} =  \frac{\partial H_y}{\partial x}-\frac{\partial H_x}{\partial y}$$
  * 3D (full field):
       * $$\frac{\partial H_x}{\partial t} = -\frac{\partial E_z}{\partial y}+\frac{\partial E_y}{\partial z}$$
       *  $$\frac{\partial H_y}{\partial t} = -\frac{\partial E_x}{\partial z}+\frac{\partial E_z}{\partial x}$$
       *  $$\frac{\partial H_z}{\partial t} = -\frac{\partial E_y}{\partial x}+\frac{\partial E_x}{\partial y}$$
       *  $$\frac{\partial E_x}{\partial t} = \frac{\partial H_z}{\partial y}-\frac{\partial H_y}{\partial z}$$
       *  $$\frac{\partial E_y}{\partial t} = \frac{\partial H_x}{\partial z}-\frac{\partial H_z}{\partial x}$$
       *  $$\frac{\partial E_z}{\partial t} = \frac{\partial H_y}{\partial x}-\frac{\partial H_x}{\partial y}$$
