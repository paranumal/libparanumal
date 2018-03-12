
/* Homogeneous Dirichlet boundary condition   */
#define ellipticDirichletCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB) \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
    uzB = uzM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticNeumannCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB) \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
    uzB = 0.f;   \
  }
