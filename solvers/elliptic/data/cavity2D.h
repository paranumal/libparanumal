
/* Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = occaCos(OCCA_PI*x)*occaCos(OCCA_PI*y);   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition2D(t,x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = uM;    \
    uxB = -OCCA_PI*occaSin(OCCA_PI*x)*occaCos(OCCA_PI*y);   \
    uyB = -OCCA_PI*occaCos(OCCA_PI*x)*occaSin(OCCA_PI*y);   \
  }
