/* Wall 1, Inflow 2, Outflow 3 */

/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in elliptic solves for velocity.
/*-----------------------------------------------------------------------------------------------*/

/* Homogeneous Dirichlet boundary condition   */
#define ellipticHomogeneousDirichlet2D(uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
  }

#define ellipticHomogeneousBC2D(bc,uM,uxM,uyM,uB,uxB,uyB)  \
  {                 \
    if     (bc==1) ellipticHomogeneousDirichlet2D(uM,uxM,uyM,uB,uxB,uyB) \
    else if(bc==2) ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB) \
    else if(bc==3) ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)  \
  }

//stub
#define ellipticBoundaryConditions2D(bc,t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB) \
  { \
  }

