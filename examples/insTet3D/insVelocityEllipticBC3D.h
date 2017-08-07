/* Wall 1, Inflow 2, Outflow 3 */

/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in elliptic solves for velocity.
/*-----------------------------------------------------------------------------------------------*/

/* Homogeneous Dirichlet boundary condition   */
#define ellipticHomogeneousDirichlet3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
    uzB = uzM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticHomogeneousNeumann3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
    uzB = 0.f;   \
  }

#define ellipticHomogeneousBC3D(bc,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {                 \
    if     (bc==1) ellipticHomogeneousDirichlet3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
    else if(bc==2) ellipticHomogeneousDirichlet3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
    else if(bc==3) ellipticHomogeneousNeumann3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  }

//stub
#define ellipticBoundaryConditions3D(bc,t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB) \
  { \
}
