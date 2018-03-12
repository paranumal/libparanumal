/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define ellipticBoundaryConditions3D(bc,t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {                 \
    if     (bc==1) ellipticDirichletCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB) \
    else if(bc==2) ellipticNeumannCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
    else           ellipticNeumannCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  }


/*-----------------------------------------------------------------------------------------------*/
/* Homogeneuous Boundary conditions used in ellipticAx.
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
    uB = uM;     \
    uxB = 0.f;   \
    uyB = 0.f;   \
    uzB = 0.f;   \
  }

/* Dirichlet 1, Neumann 2, Robin 3 (defaulted to Neumann for now) */
#define ellipticHomogeneousBC3D(bc,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {                 \
    if     (bc==1) ellipticHomogeneousDirichlet3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB) \
    else if(bc==2) ellipticHomogeneousNeumann3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
    else           ellipticHomogeneousNeumann3D(uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  }


