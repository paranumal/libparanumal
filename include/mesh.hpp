/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef MESH_HPP
#define MESH_HPP 1

#include "core.hpp"
#include "settings.hpp"
#include "ogs.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

class meshSettings_t: public settings_t {
public:
  meshSettings_t(MPI_Comm& _comm);
  void report();
};

class mesh_t {
public:
  platform_t& platform;
  meshSettings_t& settings;

  occa::properties props;

  MPI_Comm comm;
  int rank, size;

  int dim;
  int Nverts, Nfaces, NfaceVertices;

  // indices of vertex nodes
  int *vertexNodes;

  int elementType;

  hlong Nnodes=0; //global number of element vertices
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;

  dlong Nelements=0;       //local element count
  hlong NelementsGlobal=0; //global element count
  hlong *EToV; // element-to-vertex connectivity
  dlong *EToE; // element-to-element connectivity
  int   *EToF; // element-to-(local)face connectivity
  int   *EToP; // element-to-partition/process connectivity
  int   *EToB; // element-to-boundary condition type

  hlong *elementInfo; //type of element

  // boundary faces
  hlong NboundaryFaces=0; // number of boundary faces
  hlong *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  // MPI halo exchange info
  halo_t *halo;            // halo exchange pointer
  halo_t *ringHalo;        // ring halo exchange pointer
  dlong NinternalElements=0; // number of elements that can update without halo exchange
  dlong NhaloElements=0;     // number of elements that cannot update without halo exchange
  dlong  totalHaloPairs=0;   // number of elements to be received in halo exchange
  dlong  totalRingElements=0;// number of elements to be received in ring halo exchange
  dlong *internalElementIds;  // list of elements that can update without halo exchange
  dlong *haloElementIds;      // list of elements to be sent in halo exchange
  occa::memory o_internalElementIds;  // list of elements that can update without halo exchange
  occa::memory o_haloElementIds;      // list of elements to be sent in halo exchange

  // CG gather-scatter info
  ogs_t *ogs;              //occa gs pointer
  hlong *globalIds;

  // list of elements that are needed for global gather-scatter
  dlong NglobalGatherElements=0;
  dlong *globalGatherElementList;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements=0;
  dlong *localGatherElementList;
  occa::memory o_localGatherElementList;

  // volumeGeometricFactors;
  dlong Nvgeo=0;
  dfloat *vgeo;

  // second order volume geometric factors
  dlong Nggeo=0;
  dfloat *ggeo;

  // volume node info
  int N=0, Nq=0, Np=0;  // N = Polynomial order, Nq=N+1, and Np = Nodes per element
  dfloat *r, *s, *t;    // coordinates of reference nodes
  dfloat *w;            // quadrature weights (1d quadrature for tensor prod elements)
  dfloat *MM, *invMM;   // reference mass matrix

  dfloat *Dr, *Ds, *Dt; // collocation differentiation matrices
  dfloat *D;            // packed collocation differentiation matrices,
                        //  or 1D derivative for quads and hexes

  dfloat *Srr,*Srs, *Srt; //element stiffness matrices
  dfloat *Sss,*Sst, *Stt;
  dfloat *S;              // packed element stiffness matrices

  dfloat *x, *y, *z;    // coordinates of physical nodes

  /* GeoData for affine mapped elements */
  /* NC: disabling until we re-add treatment of affine elements
  dfloat *EXYZ;  // element vertices for reconstructing geofacs
  dfloat *gllzw; // GLL nodes and weights
  dfloat *ggeoNoJW;
  occa::memory o_EXYZ;
  occa::memory o_gllzw;
  occa::memory o_ggeoNoJW;
  */

  // face node info
  int Nfp=0;         // number of nodes per face
  int *faceNodes;    // list of element reference interpolation nodes on element faces
  dlong *vmapM;      // list of volume nodes that are face nodes
  dlong *vmapP;      // list of volume nodes that are paired with face nodes
  dlong *mapP;       // list of surface nodes that are paired with -ve surface  nodes
  int *faceVertices; // list of mesh vertices on each face

  dfloat *LIFT; // lift matrix
  dfloat *sM;   // surface mass MM*LIFT

  dlong   Nsgeo=0;
  dfloat *sgeo;

  // cubature
  int cubN=0, cubNp=0, cubNfp=0, cubNq=0;
  dfloat *cubr, *cubs, *cubt, *cubw;    // coordinates and weights of reference cubature nodes
  dfloat *cubx, *cuby, *cubz;           // coordinates of physical cubature nodes
  dfloat *cubInterp;                    // interpolate from W&B to cubature nodes
  dfloat *cubProject;                   // projection matrix from cubature nodes to W&B nodes
  dfloat *cubD;                         // packed differentiation matrices
  dfloat *cubPDT;                       // packed weak differentiation matrices
  dfloat *cubPDrT, *cubPDsT, *cubPDtT;  // weak differentiation matrices

  dfloat *cubvgeo;  //volume geometric data at cubature points
  dfloat *cubsgeo;  //surface geometric data at cubature points
  dfloat *cubggeo;  //second type volume geometric data at cubature points

  // surface integration node info
  int    intNfp=0;    // number of integration nodes on each face
  dfloat *intr, *ints, *intw;
  dfloat *intInterp; // interp from surface node to integration nodes
  dfloat *intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty, *intz; // coordinates of suface integration nodes

  //pml lists
  dlong NnonPmlElements=0;
  dlong NpmlElements=0;

  dlong *pmlElements;
  dlong *nonPmlElements;
  dlong *pmlIds;

  //multirate lists
  int mrNlevels=0;
  int *mrLevel;
  dlong *mrNelements, *mrInterfaceNelements;
  dlong **mrElements, **mrInterfaceElements;

  //multirate pml lists
  dlong *mrNnonPmlElements, *mrNpmlElements;
  dlong **mrPmlElements, **mrNonPmlElements;
  dlong **mrPmlIds;

  // plotting info for generating field vtu
  int    plotNverts=0;    // number of vertices for each plot element
  int    plotN=0;         // degree of plot interpolation
  int    plotNq=0;        // plotNq = plotN+1
  int    plotNp=0;        // number of plot nodes per element
  int    plotNelements=0; // number of "plot elements" per element
  int    *plotEToV;       // triangulation of plot nodes
  dfloat *plotR, *plotS, *plotT; // coordinates of plot nodes in reference element
  dfloat *plotInterp;    // reference to plot node interpolation matrix

  //SEMFEM data
  int NpFEM=0, NelFEM=0;
  int *FEMEToV;
  dfloat *rFEM, *sFEM, *tFEM;
  dfloat *SEMFEMInterp;

  // occa stuff
  occa::memory o_SEMFEMInterp;
  occa::memory o_SEMFEMAnterp;

  occa::memory o_MM;  // Mass matrix
  occa::memory o_D;   // packed differentiation matricies (contains the transpose 1d D matrix for quads/hexes)
  occa::memory o_S;   // packed stiffness matricies
  occa::memory o_LIFT;// Surface lift matrix
  occa::memory o_sM;  // Surface mass

  // volume, surface, and second order geometric factors
  occa::memory o_vgeo, o_sgeo, o_ggeo;

  //face node mappings
  occa::memory o_vmapM, o_vmapP, o_mapP;

  //element boundary mappings
  occa::memory o_EToB;

  //physical coordinates
  occa::memory o_x, o_y, o_z;

  // cubature
  occa::memory o_cubInterp, o_cubProject; //cubature interpolationm and projection
  occa::memory o_cubPDT, o_cubD;          // weak cubature derivatives, and cubature derivatives
  occa::memory o_intLIFT, o_intInterp;

  //physical cubature coordinates
  occa::memory o_cubx, o_cuby, o_cubz;

  //physical surface cubature coordinates
  occa::memory o_intx, o_inty, o_intz;

  // volume, surface, and second order geometric factors at cubature points
  occa::memory o_cubvgeo, o_cubsgeo, o_cubggeo;

  //pml lists
  occa::memory o_pmlElements;
  occa::memory o_nonPmlElements;
  occa::memory o_pmlIds;

  //multirate lists
  occa::memory o_mrLevel;
  occa::memory o_mrNelements, o_mrInterfaceNelements;
  occa::memory *o_mrElements, *o_mrInterfaceElements;

  //multirate pml lists
  occa::memory *o_mrPmlElements, *o_mrNonPmlElements;
  occa::memory *o_mrPmlIds;

  occa::kernel MassMatrixKernel;

  mesh_t() = delete;
  mesh_t(platform_t& _platform, meshSettings_t& _settings,
         MPI_Comm _comm);

  virtual ~mesh_t();

  // generic mesh setup
  static mesh_t& Setup(platform_t& _platform, meshSettings_t& _settings,
                       MPI_Comm _comm);

  // box mesh
  virtual void SetupBox() = 0;

  // pml box mesh
  virtual void SetupPmlBox() = 0;

  // mesh reader
  virtual void ParallelReader(const char *fileName) = 0;

  // repartition elements in parallel
  virtual void GeometricPartition() = 0;

  /* build parallel face connectivity */
  void ParallelConnect();
  void Connect();

  // build element-boundary connectivity
  void ConnectBoundary();

  virtual void ReferenceNodes(int N) = 0;

  /* compute x,y,z coordinates of each node */
  virtual void PhysicalNodes() = 0;

  // compute geometric factors for local to physical map
  virtual void GeometricFactors() = 0;

  virtual void SurfaceGeometricFactors() = 0;

  // serial face-node to face-node connection
  virtual void ConnectFaceNodes() = 0;

  // setup halo region
  void HaloSetup();

  // setup trace halo
  void HaloRingSetup();

  // setup trace halo
  halo_t* HaloTraceSetup(int Nfields);

  /* build global connectivity in parallel */
  void ParallelConnectNodes();

  /* build global gather scatter ops */
  void ParallelGatherScatterSetup();

  //Setup PML elements
  void PmlSetup();
  void MultiRatePmlSetup();

  //Multirate partitioning
  void MultiRateSetup(dfloat *EToDT);

  // Multirate trace halo
  halo_t** MultiRateHaloTraceSetup(int Nfields);

  virtual void OccaSetup();

  virtual void CubatureSetup()=0;

  virtual void CubatureNodes()=0;

  // print out parallel partition i
  void PrintPartitionStatistics();

  virtual dfloat ElementCharacteristicLength(dlong e) = 0;

  virtual dfloat MinCharacteristicLength() = 0;

  virtual void PlotInterp(const dfloat* q, dfloat* Iq, dfloat* scratch=nullptr)=0;

  void RecursiveSpectralBisectionPartition();

  void MassMatrixApply(occa::memory& o_q, occa::memory& o_Mq);
  virtual void MassMatrixKernelSetup(int Nfields)=0;

  //create a new mesh object with the same geometry, but different degree
  mesh_t& SetupNewDegree(int Nf);

  mesh_t* SetupRingPatch();

  mesh_t* SetupSEMFEM(hlong **globalIds, int *Nfp, int **faceNodes);

  void DegreeRaiseMatrix1D(int Nc, int Nf, dfloat *P);
  void DegreeRaiseMatrixTri2D(int Nc, int Nf, dfloat *P);
  void DegreeRaiseMatrixTet3D(int Nc, int Nf, dfloat *P);

  /***************************************************************************/
  // Basic codes for generating nodes, polynomials, matrices, etc.

public:
  //1D
  static void Nodes1D(int N, dfloat *r);
  static void EquispacedNodes1D(int _N, dfloat *_r);
  static void OrthonormalBasis1D(dfloat a, int i, dfloat *P);
  static void GradOrthonormalBasis1D(dfloat a, int i, dfloat *Pr);
  static void Vandermonde1D(int N, int Npoints, dfloat *r, dfloat *V);
  static void GradVandermonde1D(int N, int Npoints, dfloat *r, dfloat *Vr);

  static void MassMatrix1D(int _Np, dfloat *V, dfloat *MM);
  static void Dmatrix1D(int _N, int NpointsIn, dfloat *_rIn, int NpointsOut, dfloat *_rOut, dfloat *_Dr);
  static void InterpolationMatrix1D(int _N,int NpointsIn, dfloat *rIn, int NpointsOut, dfloat *rOut, dfloat *I);
  static void CubatureWeakDmatrix1D(int _Nq, int _cubNq, dfloat *_cubProject, dfloat *_cubD, dfloat *_cubPDT);

  //Jacobi polynomial evaluation
  static dfloat JacobiP(dfloat a, dfloat alpha, dfloat beta, int N);
  static dfloat GradJacobiP(dfloat a, dfloat alpha, dfloat beta, int N);

  //Gauss-Legendre-Lobatto quadrature nodes
  static void JacobiGLL(int N, dfloat *x, dfloat *w=NULL);

  //Nth order Gauss-Jacobi quadrature nodes and weights
  static void JacobiGQ(dfloat alpha, dfloat beta, int N, dfloat *x, dfloat *w);

  //Tris
  static void NodesTri2D(int _N, dfloat *_r, dfloat *_s);
  static void FaceNodesTri2D(int _N, dfloat *_r, dfloat *_s, int *_faceNodes);
  static void VertexNodesTri2D(int _N, dfloat *_r, dfloat *_s, int *_vertexNodes);
  static void EquispacedNodesTri2D(int _N, dfloat *_r, dfloat *_s);
  static void EquispacedEToVTri2D(int _N, int *_EToV);
  static void SEMFEMNodesTri2D(int _N, int *_Np, dfloat **_r, dfloat **_s);
  static void SEMFEMEToVTri2D(int _N, int *_NelFEM, int **_EToV);
  static void OrthonormalBasisTri2D(dfloat a, dfloat b, int i, int j, dfloat *P);
  static void GradOrthonormalBasisTri2D(dfloat a, dfloat b, int i, int j, dfloat *Pr, dfloat *Ps);
  static void VandermondeTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *V);
  static void GradVandermondeTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *Vr, dfloat *Vs);
  static void MassMatrixTri2D(int _Np, dfloat *V, dfloat *_MM);
  static void invMassMatrixTri2D(int _Np, dfloat *V, dfloat *_invMM);
  static void DmatrixTri2D(int _N, int Npoints, dfloat *_r, dfloat *_s,
                           dfloat *_Dr, dfloat *_Ds);
  static void LIFTmatrixTri2D(int _N, int *_faceNodes,
                              dfloat *_r, dfloat *_s, dfloat *_LIFT);
  static void SurfaceMassMatrixTri2D(int _N, dfloat *_MM, dfloat *_LIFT, dfloat *_sM);
  static void SmatrixTri2D(int _N, dfloat *_Dr, dfloat *_Ds, dfloat *_MM,
                           dfloat *_Srr, dfloat *_Srs, dfloat *_Sss);
  static void InterpolationMatrixTri2D(int _N,
                                       int NpointsIn, dfloat *rIn, dfloat *sIn,
                                       int NpointsOut, dfloat *rOut, dfloat *sOut,
                                       dfloat *I);
  static void CubatureNodesTri2D(int cubTriN, int*cubNp, dfloat **cubTrir, dfloat **cubTris, dfloat **cubTriw);
  static void CubaturePmatrixTri2D(int _N, int _Np, dfloat *_r, dfloat *_s,
                                   int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubProject);
  static void CubatureWeakDmatricesTri2D(int _N, int _Np, dfloat *_r, dfloat *_s,
                                         int _cubNp, dfloat *_cubr, dfloat *_cubs,
                                         dfloat *_cubPDrT, dfloat *_cubPDsT);
  static void CubatureSurfaceMatricesTri2D(int _N, int _Np, dfloat *_r, dfloat *_s, int *_faceNodes,
                                           int _intNfp, dfloat *_intr, dfloat *_intw,
                                           dfloat *_intInterp, dfloat *_intLIFT);
  static void SEMFEMInterpMatrixTri2D(int _N,
                                      int _Np, dfloat *_r, dfloat *_s,
                                      int _NpFEM, dfloat *rFEM, dfloat *sFEM,
                                      dfloat *I);

  static void Warpfactor(int _N, int Npoints, dfloat *r, dfloat *w);
  static void WarpBlendTransformTri2D(int _N, int _Npoints, dfloat *_r, dfloat *_s, dfloat alphaIn=-1);


  //Quads
  static void NodesQuad2D(int _N, dfloat *_r, dfloat *_s);
  static void FaceNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_faceNodes);
  static void VertexNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_vertexNodes);
  static void EquispacedNodesQuad2D(int _N, dfloat *_r, dfloat *_s);
  static void EquispacedEToVQuad2D(int _N, int *_EToV);
  static void SEMFEMEToVQuad2D(int _N, int *_EToV);
  static void OrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *P);
  static void GradOrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *Pr, dfloat *Ps);
  static void VandermondeQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *V);
  static void GradVandermondeQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *Vr, dfloat *Vs);
  static void MassMatrixQuad2D(int _Np, dfloat *V, dfloat *_MM);
  static void LumpedMassMatrixQuad2D(int _N, dfloat *_gllw, dfloat *_MM);
  static void invLumpedMassMatrixQuad2D(int _N, dfloat *_gllw, dfloat *_invMM);
  static void DmatrixQuad2D(int _N, int Npoints, dfloat *_r, dfloat *_s,
                                          dfloat *_Dr, dfloat *_Ds);
  static void InterpolationMatrixQuad2D(int _N,
                                        int NpointsIn, dfloat *rIn, dfloat *sIn,
                                        int NpointsOut, dfloat *rOut, dfloat *sOut,
                                        dfloat *I);

  //Tets
  static void NodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  static void FaceNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes);
  static void VertexNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_vertexNodes);
  static void EquispacedNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  static void EquispacedEToVTet3D(int _N, int *_EToV);
  static void SEMFEMEToVTet3D(int _N, int *_EToV);
  static void OrthonormalBasisTet3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P);
  static void GradOrthonormalBasisTet3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *Pr, dfloat *Ps, dfloat *Pt);
  static void VandermondeTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *V);
  static void GradVandermondeTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *Vr, dfloat *Vs, dfloat *Vt);
  static void MassMatrixTet3D(int _Np, dfloat *V, dfloat *_MM);
  static void invMassMatrixTet3D(int _Np, dfloat *V, dfloat *_invMM);
  static void DmatrixTet3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                           dfloat *_Dr, dfloat *_Ds, dfloat *_Dt);
  static void LIFTmatrixTet3D(int _N, int *_faceNodes,
                              dfloat *_r, dfloat *_s, dfloat *_t, dfloat *_LIFT);
  static void SurfaceMassMatrixTet3D(int _N, dfloat *_MM, dfloat *_LIFT, dfloat *_sM);
  static void SmatrixTet3D(int _N, dfloat *_Dr, dfloat *_Ds, dfloat *_Dt, dfloat *_MM,
                           dfloat *_Srr, dfloat *_Srs, dfloat *_Srt,
                           dfloat *_Sss, dfloat *_Sst, dfloat *_Stt);
  static void InterpolationMatrixTet3D(int _N,
                                       int NpointsIn, dfloat *rIn, dfloat *sIn, dfloat *tIn,
                                       int NpointsOut, dfloat *rOut, dfloat *sOut, dfloat *tOut,
                                       dfloat *I);
  static void CubatureNodesTet3D(int cubN, int*cubNp, dfloat **cubr, dfloat **cubs, dfloat **cubt, dfloat **cubw);
  static void CubaturePmatrixTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                   int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubt,
                                   dfloat *_cubProject);
  static void CubatureWeakDmatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                         int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubt,
                                         dfloat *_cubPDrT, dfloat *_cubPDsT, dfloat *_cubPDtT);
  static void CubatureSurfaceMatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes,
                                           int _intNfp, dfloat *_intr, dfloat *_ints, dfloat *_intw,
                                           dfloat *_intInterp, dfloat *_intLIFT);
  static void SEMFEMInterpMatrixTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                      int _NpFEM, dfloat *_rFEM, dfloat *_sFEM, dfloat *_tFEM,
                                      dfloat *I);
  static void WarpShiftFace3D(int _N, int Npoints, dfloat alpha,
                              dfloat *L1, dfloat *L2, dfloat *L3,
                              dfloat *w1, dfloat *w2);
  static void WarpBlendTransformTet3D(int _N, int _Npoints, dfloat *_r, dfloat *_s, dfloat *_t, dfloat alphaIn=-1);


  //Hexs
  static void NodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  static void FaceNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t,  int *_faceNodes);
  static void VertexNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_vertexNodes);
  static void EquispacedNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  static void EquispacedEToVHex3D(int _N, int *_EToV);
  static void SEMFEMEToVHex3D(int _N, int *_EToV);
  static void OrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P);
  static void GradOrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *Pr, dfloat *Ps, dfloat *Pt);
  static void VandermondeHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *V);
  static void GradVandermondeHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t,
                                   dfloat *Vr, dfloat *Vs, dfloat *Vt);
  static void MassMatrixHex3D(int _Np, dfloat *V, dfloat *_MM);
  static void LumpedMassMatrixHex3D(int _N, dfloat *_gllw, dfloat *_MM);
  static void invLumpedMassMatrixHex3D(int _N, dfloat *_gllw, dfloat *_invMM);
  static void DmatrixHex3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                           dfloat *_Dr, dfloat *_Ds, dfloat *_Dt);
  static void InterpolationMatrixHex3D(int _N,
                                       int NpointsIn, dfloat *rIn, dfloat *sIn, dfloat *tIn,
                                       int NpointsOut, dfloat *rOut, dfloat *sOut, dfloat *tOut,
                                       dfloat *I);
};

#endif

