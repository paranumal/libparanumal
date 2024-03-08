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
