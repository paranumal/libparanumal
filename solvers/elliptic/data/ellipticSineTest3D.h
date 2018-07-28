
/* Homogeneous Dirichlet boundary condition   */
#define ellipticDirichletCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    uB  = cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);	\
    uxB = uxM;	\
    uyB = uyM;   \
    uzB = uzM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticNeumannCondition3D(t,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    uB  = uM;    \
    uxB = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);	\
    uyB = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);		\
    uzB = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);		\
  }

