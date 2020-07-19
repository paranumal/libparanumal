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

#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#include "types.h"
#include "utils.hpp"
#include "ogs.hpp"
#include "settings.hpp"

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
  occa::device& device;
  MPI_Comm& comm;
  meshSettings_t& settings;
  occa::properties& props;

  int rank, size; // MPI rank and size (process count)

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



  mesh_t() = delete;
  mesh_t(occa::device& device, MPI_Comm& comm,
         meshSettings_t& settings, occa::properties& props);

  virtual ~mesh_t();

  // generic mesh setup
  static mesh_t& Setup(occa::device& device, MPI_Comm& comm,
                       meshSettings_t& settings, occa::properties& props);

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

  void RecursiveSpectralBisectionPartition();

  //create a new mesh object with the same geometry, but different degree
  mesh_t& SetupNewDegree(int Nf);

  mesh_t* SetupRingPatch();

  mesh_t* SetupSEMFEM(hlong **globalIds, int *Nfp, int **faceNodes);

  void DegreeRaiseMatrix1D(int Nc, int Nf, dfloat *P);
  void DegreeRaiseMatrixTri2D(int Nc, int Nf, dfloat *P);
  void DegreeRaiseMatrixTet3D(int Nc, int Nf, dfloat *P);
protected:
  //1D
  void Nodes1D(int N, dfloat *r);
  void EquispacedNodes1D(int _N, dfloat *_r);
  void OrthonormalBasis1D(dfloat a, int i, dfloat *P);
  void GradOrthonormalBasis1D(dfloat a, int i, dfloat *Pr);
  void Vandermonde1D(int N, int Npoints, dfloat *r, dfloat *V);
  void GradVandermonde1D(int N, int Npoints, dfloat *r, dfloat *Vr);

  void MassMatrix1D(int _Np, dfloat *V, dfloat *MM);
  void Dmatrix1D(int _N, int NpointsIn, dfloat *_rIn, int NpointsOut, dfloat *_rOut, dfloat *_Dr);
  void InterpolationMatrix1D(int _N,int NpointsIn, dfloat *rIn, int NpointsOut, dfloat *rOut, dfloat *I);
  void CubatureWeakDmatrix1D(int _Nq, int _cubNq, dfloat *_cubProject, dfloat *_cubD, dfloat *_cubPDT);

  //Jacobi polynomial evaluation
  dfloat JacobiP(dfloat a, dfloat alpha, dfloat beta, int N);
  dfloat GradJacobiP(dfloat a, dfloat alpha, dfloat beta, int N);

  //Gauss-Legendre-Lobatto quadrature nodes
  void JacobiGLL(int N, dfloat *x, dfloat *w=NULL);

  //Nth order Gauss-Jacobi quadrature nodes and weights
  void JacobiGQ(dfloat alpha, dfloat beta, int N, dfloat *x, dfloat *w);

  //Tris
  void NodesTri2D(int _N, dfloat *_r, dfloat *_s);
  void FaceNodesTri2D(int _N, dfloat *_r, dfloat *_s, int *_faceNodes);
  void VertexNodesTri2D(int _N, dfloat *_r, dfloat *_s, int *_vertexNodes);
  void EquispacedNodesTri2D(int _N, dfloat *_r, dfloat *_s);
  void EquispacedEToVTri2D(int _N, int *_EToV);
  void SEMFEMNodesTri2D(int _N, int *_Np, dfloat **_r, dfloat **_s);
  void SEMFEMEToVTri2D(int _N, int *_NelFEM, int **_EToV);
  void OrthonormalBasisTri2D(dfloat a, dfloat b, int i, int j, dfloat *P);
  void GradOrthonormalBasisTri2D(dfloat a, dfloat b, int i, int j, dfloat *Pr, dfloat *Ps);
  void VandermondeTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *V);
  void GradVandermondeTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *Vr, dfloat *Vs);
  void MassMatrixTri2D(int _Np, dfloat *V, dfloat *_MM);
  void invMassMatrixTri2D(int _Np, dfloat *V, dfloat *_invMM);
  void DmatrixTri2D(int _N, int Npoints, dfloat *_r, dfloat *_s,
                                          dfloat *_Dr, dfloat *_Ds);
  void LIFTmatrixTri2D(int _N, int *_faceNodes,
                             dfloat *_r, dfloat *_s, dfloat *_LIFT);
  void SurfaceMassMatrixTri2D(int _N, dfloat *_MM, dfloat *_LIFT, dfloat *_sM);
  void SmatrixTri2D(int _N, dfloat *_Dr, dfloat *_Ds, dfloat *_MM,
                          dfloat *_Srr, dfloat *_Srs, dfloat *_Sss);
  void InterpolationMatrixTri2D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut,
                               dfloat *I);
  void CubatureNodesTri2D(int cubTriN, int*cubNp, dfloat **cubTrir, dfloat **cubTris, dfloat **cubTriw);
  void CubaturePmatrixTri2D(int _N, int _Np, dfloat *_r, dfloat *_s,
                            int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubProject);
  void CubatureWeakDmatricesTri2D(int _N, int _Np, dfloat *_r, dfloat *_s,
                                  int _cubNp, dfloat *_cubr, dfloat *_cubs,
                                  dfloat *_cubPDrT, dfloat *_cubPDsT);
  void CubatureSurfaceMatricesTri2D(int _N, int _Np, dfloat *_r, dfloat *_s, int *_faceNodes,
                                    int _intNfp, dfloat *_intr, dfloat *_intw,
                                    dfloat *_intInterp, dfloat *_intLIFT);
  void SEMFEMInterpMatrixTri2D(int _N,
                                int _Np, dfloat *_r, dfloat *_s,
                                int _NpFEM, dfloat *rFEM, dfloat *sFEM,
                                dfloat *I);

  void Warpfactor(int _N, int Npoints, dfloat *r, dfloat *w);
  void WarpBlendTransformTri2D(int _N, int _Npoints, dfloat *_r, dfloat *_s, dfloat alphaIn=-1);


  //Quads
  void NodesQuad2D(int _N, dfloat *_r, dfloat *_s);
  void FaceNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_faceNodes);
  void VertexNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_vertexNodes);
  void EquispacedNodesQuad2D(int _N, dfloat *_r, dfloat *_s);
  void EquispacedEToVQuad2D(int _N, int *_EToV);
  void SEMFEMEToVQuad2D(int _N, int *_EToV);
  void OrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *P);
  void GradOrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *Pr, dfloat *Ps);
  void VandermondeQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *V);
  void GradVandermondeQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *Vr, dfloat *Vs);
  void MassMatrixQuad2D(int _Np, dfloat *V, dfloat *_MM);
  void LumpedMassMatrixQuad2D(int _N, dfloat *_gllw, dfloat *_MM);
  void invLumpedMassMatrixQuad2D(int _N, dfloat *_gllw, dfloat *_invMM);
  void DmatrixQuad2D(int _N, int Npoints, dfloat *_r, dfloat *_s,
                                          dfloat *_Dr, dfloat *_Ds);
  void InterpolationMatrixQuad2D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut,
                               dfloat *I);

  //Tets
  void NodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  void FaceNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes);
  void VertexNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_vertexNodes);
  void EquispacedNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  void EquispacedEToVTet3D(int _N, int *_EToV);
  void SEMFEMEToVTet3D(int _N, int *_EToV);
  void OrthonormalBasisTet3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P);
  void GradOrthonormalBasisTet3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *Pr, dfloat *Ps, dfloat *Pt);
  void VandermondeTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *V);
  void GradVandermondeTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *Vr, dfloat *Vs, dfloat *Vt);
  void MassMatrixTet3D(int _Np, dfloat *V, dfloat *_MM);
  void invMassMatrixTet3D(int _Np, dfloat *V, dfloat *_invMM);
  void DmatrixTet3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                                          dfloat *_Dr, dfloat *_Ds, dfloat *_Dt);
  void LIFTmatrixTet3D(int _N, int *_faceNodes,
                             dfloat *_r, dfloat *_s, dfloat *_t, dfloat *_LIFT);
  void SurfaceMassMatrixTet3D(int _N, dfloat *_MM, dfloat *_LIFT, dfloat *_sM);
  void SmatrixTet3D(int _N, dfloat *_Dr, dfloat *_Ds, dfloat *_Dt, dfloat *_MM,
                          dfloat *_Srr, dfloat *_Srs, dfloat *_Srt,
                          dfloat *_Sss, dfloat *_Sst, dfloat *_Stt);
  void InterpolationMatrixTet3D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn, dfloat *tIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut, dfloat *tOut,
                               dfloat *I);
  void CubatureNodesTet3D(int cubN, int*cubNp, dfloat **cubr, dfloat **cubs, dfloat **cubt, dfloat **cubw);
  void CubaturePmatrixTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                          int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubt,
                          dfloat *_cubProject);
  void CubatureWeakDmatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubt,
                                dfloat *_cubPDrT, dfloat *_cubPDsT, dfloat *_cubPDtT);
  void CubatureSurfaceMatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes,
                                    int _intNfp, dfloat *_intr, dfloat *_ints, dfloat *_intw,
                                    dfloat *_intInterp, dfloat *_intLIFT);
  void SEMFEMInterpMatrixTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                    int _NpFEM, dfloat *_rFEM, dfloat *_sFEM, dfloat *_tFEM,
                                    dfloat *I);
  void WarpShiftFace3D(int _N, int Npoints, dfloat alpha,
                             dfloat *L1, dfloat *L2, dfloat *L3,
                             dfloat *w1, dfloat *w2);
  void WarpBlendTransformTet3D(int _N, int _Npoints, dfloat *_r, dfloat *_s, dfloat *_t, dfloat alphaIn=-1);


  //Hexs
  void NodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  void FaceNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t,  int *_faceNodes);
  void VertexNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_vertexNodes);
  void EquispacedNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t);
  void EquispacedEToVHex3D(int _N, int *_EToV);
  void SEMFEMEToVHex3D(int _N, int *_EToV);
  void OrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P);
  void GradOrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *Pr, dfloat *Ps, dfloat *Pt);
  void VandermondeHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat *V);
  void GradVandermondeHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t,
                                                dfloat *Vr, dfloat *Vs, dfloat *Vt);
  void MassMatrixHex3D(int _Np, dfloat *V, dfloat *_MM);
  void LumpedMassMatrixHex3D(int _N, dfloat *_gllw, dfloat *_MM);
  void invLumpedMassMatrixHex3D(int _N, dfloat *_gllw, dfloat *_invMM);
  void DmatrixHex3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                                          dfloat *_Dr, dfloat *_Ds, dfloat *_Dt);
  void InterpolationMatrixHex3D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn, dfloat *tIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut, dfloat *tOut,
                               dfloat *I);
};

#endif

