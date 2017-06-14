/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define ellipticBoundaryConditions2D(bc,t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {                 \
    if     (bc==1) ellipticDiricletCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB) \
    else if(bc==2) ellipticNeumannCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
    else           ellipticNeumannCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  }


/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in ellipticAx.
/*-----------------------------------------------------------------------------------------------*/

/* Homogeneous Dirichlet boundary condition   */
#define ellipticHomogeneousDiriclet2D(uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB = uM;     \
    uxB = 0.f;   \
    uyB = 0.f;   \
  }

/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define ellipticHomogeneousBC2D(bc,uM,uxM,uyM,uB,uxB,uyB)  \
  {                 \
    if     (bc==1) ellipticHomogeneousDiriclet2D(uM,uxM,uyM,uB,uxB,uyB) \
    else if(bc==2) ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)  \
    else           ellipticHomogeneousNeumann2D(uM,uxM,uyM,uB,uxB,uyB)  \
  }

