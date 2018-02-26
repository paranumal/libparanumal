
/* Homogeneous Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticNeumannCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
  }
