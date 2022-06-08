/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
#include "platform.hpp"
#include "settings.hpp"
#include "ogs.hpp"

namespace libp {

class meshSettings_t: public settings_t {
public:
  meshSettings_t() = default;
  meshSettings_t(comm_t _comm);
  void report();
};

namespace Mesh {
  /*Element types*/
  enum ElementType {
    TRIANGLES     =3,
    QUADRILATERALS=4,
    TETRAHEDRA    =6,
    HEXAHEDRA     =12
  };
} //namespace Mesh

class mesh_t {
 public:
  platform_t platform;
  meshSettings_t settings;
  properties_t props;

  comm_t comm;
  int rank, size;

  /*************************/
  /* Element Data          */
  /*************************/
  int dim;
  int Nverts, Nfaces, NfaceVertices;
  Mesh::ElementType elementType;

  // indices of vertex nodes
  memory<int> vertexNodes;

  hlong Nnodes=0; //global number of element vertices
  memory<dfloat> EX; // coordinates of vertices for each element
  memory<dfloat> EY;
  memory<dfloat> EZ;

  dlong Nelements=0;       //local element count
  hlong NelementsGlobal=0; //global element count
  memory<hlong> EToV;      // element-to-vertex connectivity
  memory<hlong> EToE;      // element-to-element connectivity
  memory<int>   EToF;      // element-to-(local)face connectivity
  memory<int>   EToP;      // element-to-partition/process connectivity
  memory<int>   EToB;      // element-to-boundary condition type
  deviceMemory<int> o_EToB;

  memory<int>   mapB;      // node-to-boundary condition type
  deviceMemory<int> o_mapB;

  memory<hlong> elementInfo; //type of element

  memory<dlong> VmapM;  // list of vertices on each face
  memory<dlong> VmapP;  // list of vertices that are paired with face vertices

  // boundary faces
  hlong NboundaryFaces=0; // number of boundary faces
  memory<hlong> boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  /*************************/
  /* FEM Space             */
  /*************************/
  int N=0, Np=0;             // N = Polynomial order and Np = Nodes per element
  memory<dfloat> r, s, t;    // coordinates of local nodes

  int Nq=0;                 // N = Polynomial order, Nq=N+1
  memory<dfloat> gllz;      // 1D GLL quadrature nodes
  memory<dfloat> gllw;      // 1D GLL quadrature weights

  // face node info
  int Nfp=0;                // number of nodes per face
  memory<int> faceNodes;    // list of element reference interpolation nodes on element faces
  memory<int> faceVertices; // list of mesh vertices on each face

  /*************************/
  /* FEM Operators         */
  /*************************/
  memory<dfloat> Dr, Ds, Dt;    // collocation differentiation matrices
  memory<dfloat> D;
  deviceMemory<dfloat> o_D;
  memory<dfloat> MM, invMM;     // reference mass matrix
  deviceMemory<dfloat> o_MM;
  memory<dfloat> LIFT;          // lift matrix
  deviceMemory<dfloat> o_LIFT;
  memory<dfloat> sM;            // surface mass (MM*LIFT)^T
  deviceMemory<dfloat> o_sM;
  memory<dfloat> Srr, Srs, Srt; //element stiffness matrices
  memory<dfloat> Ssr, Sss, Sst;
  memory<dfloat> Str, Sts, Stt;
  memory<dfloat> S;
  deviceMemory<dfloat> o_S;

  /*************************/
  /* Cubature              */
  /*************************/
  // cubature
  int cubN=0, cubNp=0, cubNfp=0, cubNq=0;
  memory<dfloat> cubr, cubs, cubt, cubw; // coordinates and weights of local cubature nodes

  memory<dfloat> cubInterp;    // interpolate from W&B to cubature nodes
  deviceMemory<dfloat> o_cubInterp;
  memory<dfloat> cubProject;   // projection matrix from cubature nodes to W&B nodes
  deviceMemory<dfloat> o_cubProject;
  memory<dfloat> cubD;         // 1D differentiation matrix
  deviceMemory<dfloat> o_cubD;
  memory<dfloat> cubPDrT, cubPDsT, cubPDtT;  // weak differentiation matrices
  memory<dfloat> cubPDT;                     // packed weak differentiation matrices
  deviceMemory<dfloat> o_cubPDT;

  // surface integration node info
  int intNfp=0;    // number of integration nodes on each face
  memory<dfloat> intr, ints, intw;
  memory<dfloat> intInterp; // interp from surface node to integration nodes
  deviceMemory<dfloat> o_intInterp;
  memory<dfloat> intLIFT;   // lift from surface integration nodes to W&B volume nodes
  deviceMemory<dfloat> o_intLIFT;

  /*************************/
  /* Plotting              */
  /*************************/
  // ploting info for generating field vtu
  int plotN=0;
  int plotNq=0;
  int plotNp=0;
  int plotNverts;    // number of vertices for each plot element
  int plotNelements; // number of "plot elements" per element
  memory<int>   plotEToV;             // triangulation of plot nodes
  memory<dfloat> plotR, plotS, plotT; // coordinates of plot nodes in reference element
  memory<dfloat> plotInterp;          // reference to plot node interpolation matrix

  /*************************/
  /* Physical Space        */
  /*************************/
  // volume node info
  memory<dfloat> x, y, z;    // coordinates of physical nodes
  deviceMemory<dfloat> o_x, o_y, o_z;    // coordinates of physical nodes

  memory<dlong> vmapM;      // list of volume nodes that are face nodes
  deviceMemory<dlong> o_vmapM;
  memory<dlong> vmapP;      // list of volume nodes that are paired with face nodes
  deviceMemory<dlong> o_vmapP;
  memory<dlong> mapP;       // list of surface nodes that are paired with -ve surface  nodes
  deviceMemory<dlong> o_mapP;

  // Jacobian
  memory<dfloat> wJ;
  deviceMemory<dfloat> o_wJ;
  // volumeGeometricFactors;
  dlong Nvgeo;
  memory<dfloat> vgeo;
  deviceMemory<dfloat> o_vgeo;
  // surfaceGeometricFactors;
  dlong   Nsgeo;
  memory<dfloat> sgeo;
  deviceMemory<dfloat> o_sgeo;
  // second order volume geometric factors
  dlong Nggeo;
  memory<dfloat> ggeo;
  deviceMemory<dfloat> o_ggeo;

  memory<dfloat> cubx, cuby, cubz; // coordinates of physical nodes
  deviceMemory<dfloat> o_cubx, o_cuby, o_cubz;
  memory<dfloat> intx, inty, intz; // coordinates of suface integration nodes
  deviceMemory<dfloat> o_intx, o_inty, o_intz;

  memory<dfloat> cubwJ;            //Jacobian at cubature points
  deviceMemory<dfloat> o_cubwJ;
  memory<dfloat> cubvgeo;          //volume geometric data at cubature points
  deviceMemory<dfloat> o_cubvgeo;
  memory<dfloat> cubsgeo;          //surface geometric data at cubature points
  deviceMemory<dfloat> o_cubsgeo;
  memory<dfloat> cubggeo;          //second type volume geometric data at cubature points
  deviceMemory<dfloat> o_cubggeo;

  /*************************/
  /* MPI Data              */
  /*************************/
  // MPI halo exchange info
  ogs::halo_t halo;            // halo exchange pointer
  ogs::halo_t ringHalo;        // ring halo exchange pointer
  dlong NinternalElements=0; // number of elements that can update without halo exchange
  dlong NhaloElements=0;     // number of elements that cannot update without halo exchange
  dlong  totalHaloPairs=0;   // number of elements to be received in halo exchange
  dlong  totalRingElements=0;// number of elements to be received in ring halo exchange

  memory<dlong> internalElementIds;  // list of elements that can update without halo exchange
  memory<dlong> haloElementIds;      // list of elements to be sent in halo exchange
  deviceMemory<dlong> o_internalElementIds;  // list of elements that can update without halo exchange
  deviceMemory<dlong> o_haloElementIds;      // list of elements to be sent in halo exchange

  // CG gather-scatter info
  ogs::ogs_t ogs;              //occa gs pointer
  memory<hlong> globalIds;

  // list of elements that are needed for global gather-scatter
  dlong NglobalGatherElements;
  memory<dlong> globalGatherElementList;
  deviceMemory<dlong> o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements;
  memory<dlong> localGatherElementList;
  deviceMemory<dlong> o_localGatherElementList;

  /*************************/
  /* PML                   */
  /*************************/
  //pml lists
  dlong NnonPmlElements=0;
  dlong NpmlElements=0;

  memory<dlong> pmlElements;
  deviceMemory<dlong> o_pmlElements;
  memory<dlong> nonPmlElements;
  deviceMemory<dlong> o_nonPmlElements;
  memory<dlong> pmlIds;
  deviceMemory<dlong> o_pmlIds;


  /*************************/
  /* Multirate timestepping*/
  /*************************/
  //multirate lists
  int mrNlevels=0;
  memory<int> mrLevel;
  deviceMemory<int> o_mrLevel;

  memory<dlong> mrNelements, mrInterfaceNelements;
  deviceMemory<dlong> o_mrNelements, o_mrInterfaceNelements;

  memory<dlong> mrNnonPmlElements, mrNpmlElements;

  memory<memory<dlong>> mrElements, mrInterfaceElements;
  memory<deviceMemory<dlong>> o_mrElements, o_mrInterfaceElements;

  //multirate pml lists
  memory<memory<dlong>> mrPmlElements, mrNonPmlElements;
  memory<deviceMemory<dlong>> o_mrPmlElements, o_mrNonPmlElements;
  memory<memory<dlong>> mrPmlIds;
  memory<deviceMemory<dlong>> o_mrPmlIds;

  /*************************/
  /* SEMFEM                */
  /*************************/
  //SEMFEM data
  int NpFEM=0, NelFEM=0;
  memory<int> FEMEToV;
  memory<dfloat> rFEM, sFEM, tFEM;
  memory<dfloat> SEMFEMInterp;
  deviceMemory<dfloat> o_SEMFEMInterp;
  deviceMemory<dfloat> o_SEMFEMAnterp;

  kernel_t MassMatrixKernel;

  mesh_t() = default;
  mesh_t(platform_t& _platform, meshSettings_t& _settings,
         comm_t _comm) {
    Setup(_platform, _settings, _comm);
  }

  // mesh setup
  void Setup(platform_t& _platform, meshSettings_t& _settings,
             comm_t _comm);

  // setup trace halo
  void HaloRingSetup();

  // setup trace halo
  ogs::halo_t HaloTraceSetup(int Nfields);

  //Setup PML elements
  void PmlSetup();
  void MultiRatePmlSetup();

  //Multirate partitioning
  void MultiRateSetup(memory<dfloat> EToDT);

  // Multirate trace halo
  memory<ogs::halo_t> MultiRateHaloTraceSetup(int Nfields);

  // Setup cubature
  void CubatureSetup() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        CubatureSetupTri2D();
        break;
      case Mesh::QUADRILATERALS:
        CubatureSetupQuad2D();
        break;
      case Mesh::TETRAHEDRA:
        CubatureSetupTet3D();
        break;
      case Mesh::HEXAHEDRA:
        CubatureSetupHex3D();
        break;
    }
  }

  // Setup cubature physical nodes
  void CubaturePhysicalNodes() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        if (dim==2)
          CubaturePhysicalNodesTri2D();
        else
          CubaturePhysicalNodesTri3D();
        break;
      case Mesh::QUADRILATERALS:
        if (dim==2)
          CubaturePhysicalNodesQuad2D();
        else
          CubaturePhysicalNodesQuad3D();
        break;
      case Mesh::TETRAHEDRA:
        CubaturePhysicalNodesTet3D();
        break;
      case Mesh::HEXAHEDRA:
        CubaturePhysicalNodesHex3D();
        break;
    }
  }

  dfloat MinCharacteristicLength();

  void PlotInterp(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch=memory<dfloat>()) {
    switch (elementType) {
      case Mesh::TRIANGLES:
        PlotInterpTri2D(q, Iq, scratch);
        break;
      case Mesh::QUADRILATERALS:
        PlotInterpQuad2D(q, Iq, scratch);
        break;
      case Mesh::TETRAHEDRA:
        PlotInterpTet3D(q, Iq, scratch);
        break;
      case Mesh::HEXAHEDRA:
        PlotInterpHex3D(q, Iq, scratch);
        break;
    }
  }

  void MassMatrixApply(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_Mq);
  void MassMatrixKernelSetup(int Nfields) {
    switch (elementType) {
      case Mesh::TRIANGLES:
        MassMatrixKernelSetupTri2D(Nfields);
        break;
      case Mesh::QUADRILATERALS:
        MassMatrixKernelSetupQuad2D(Nfields);
        break;
      case Mesh::TETRAHEDRA:
        MassMatrixKernelSetupTet3D(Nfields);
        break;
      case Mesh::HEXAHEDRA:
        MassMatrixKernelSetupHex3D(Nfields);
        break;
    }
  }

  dfloat ElementCharacteristicLength(dlong e) {
    switch (elementType) {
      case Mesh::TRIANGLES:
        return ElementCharacteristicLengthTri2D(e);
      case Mesh::QUADRILATERALS:
        return ElementCharacteristicLengthQuad2D(e);
      case Mesh::TETRAHEDRA:
        return ElementCharacteristicLengthTet3D(e);
      case Mesh::HEXAHEDRA:
        return ElementCharacteristicLengthHex3D(e);
      default:
        return 0.0;
    }
  }

  //create a new mesh object with the same geometry, but different degree
  mesh_t SetupNewDegree(int Nf);

  mesh_t SetupRingPatch();

  mesh_t SetupSEMFEM(memory<hlong>& globalIds, memory<int>& mapB);

  int RXID, RYID, RZID;
  int SXID, SYID, SZID;
  int TXID, TYID, TZID;
  int JID, JWID, IJWID;
  int G00ID, G01ID, G02ID, G11ID, G12ID, G22ID;

  int NXID, NYID, NZID;
  int SJID, IJID, IHID, WIJID, WSJID;

 private:
  /*Set the type of mesh*/
  void SetElementType(const Mesh::ElementType eType);

  // box mesh
  void SetupBox() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        SetupBoxTri2D();
        break;
      case Mesh::QUADRILATERALS:
        SetupBoxQuad2D();
        break;
      case Mesh::TETRAHEDRA:
        SetupBoxTet3D();
        break;
      case Mesh::HEXAHEDRA:
        SetupBoxHex3D();
        break;
    }
  }
  void SetupBoxTri2D();
  void SetupBoxQuad2D();
  void SetupBoxTet3D();
  void SetupBoxHex3D();

  // pml box mesh
  void SetupPmlBox() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        SetupPmlBoxTri2D();
        break;
      case Mesh::QUADRILATERALS:
        SetupPmlBoxQuad2D();
        break;
      case Mesh::TETRAHEDRA:
        SetupPmlBoxTet3D();
        break;
      case Mesh::HEXAHEDRA:
        SetupPmlBoxHex3D();
        break;
    }
  }
  void SetupPmlBoxTri2D();
  void SetupPmlBoxQuad2D();
  void SetupPmlBoxTet3D();
  void SetupPmlBoxHex3D();

  // mesh reader
  void ReadGmsh(const std::string fileName) {
    switch (elementType) {
      case Mesh::TRIANGLES:
        if(dim==2)
          ReadGmshTri2D(fileName);
        else
          ReadGmshTri3D(fileName);
        break;
      case Mesh::QUADRILATERALS:
        if(dim==2)
          ReadGmshQuad2D(fileName);
        else
          ReadGmshQuad3D(fileName);
        break;
      case Mesh::TETRAHEDRA:
        ReadGmshTet3D(fileName);
        break;
      case Mesh::HEXAHEDRA:
        ReadGmshHex3D(fileName);
        break;
    }
  }
  void ReadGmshTri2D(const std::string fileName);
  void ReadGmshTri3D(const std::string fileName);
  void ReadGmshQuad2D(const std::string fileName);
  void ReadGmshQuad3D(const std::string fileName);
  void ReadGmshTet3D(const std::string fileName);
  void ReadGmshHex3D(const std::string fileName);

  // reference nodes and operators
  void ReferenceNodes() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        ReferenceNodesTri2D();
        break;
      case Mesh::QUADRILATERALS:
        ReferenceNodesQuad2D();
        break;
      case Mesh::TETRAHEDRA:
        ReferenceNodesTet3D();
        break;
      case Mesh::HEXAHEDRA:
        ReferenceNodesHex3D();
        break;
    }
  }
  void ReferenceNodesTri2D();
  void ReferenceNodesQuad2D();
  void ReferenceNodesTet3D();
  void ReferenceNodesHex3D();

  // repartition elements
  void Partition();

  /* build parallel face connectivity */
  void Connect();

  // build element-boundary connectivity
  void ConnectBoundary();

  // face-vertex to face-vertex connection
  void ConnectFaceVertices();

  // face-node to face-node connection
  void ConnectFaceNodes();

  // setup halo region
  void HaloSetup();

  /* build global connectivity in parallel */
  void ConnectNodes();

  /* build global gather scatter ops */
  void GatherScatterSetup();

  /* compute x,y,z coordinates of each node */
  void PhysicalNodes() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        if(dim==2)
          PhysicalNodesTri2D();
        else
          PhysicalNodesTri3D();
        break;
      case Mesh::QUADRILATERALS:
        if(dim==2)
          PhysicalNodesQuad2D();
        else
          PhysicalNodesQuad3D();
        break;
      case Mesh::TETRAHEDRA:
        PhysicalNodesTet3D();
        break;
      case Mesh::HEXAHEDRA:
        PhysicalNodesHex3D();
        break;
    }
  }
  void PhysicalNodesTri2D();
  void PhysicalNodesTri3D();
  void PhysicalNodesQuad2D();
  void PhysicalNodesQuad3D();
  void PhysicalNodesTet3D();
  void PhysicalNodesHex3D();

  // compute geometric factors for local to physical map
  void GeometricFactors() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        if(dim==2)
          GeometricFactorsTri2D();
        else
          GeometricFactorsTri3D();
        break;
      case Mesh::QUADRILATERALS:
        if(dim==2)
          GeometricFactorsQuad2D();
        else
          GeometricFactorsQuad3D();
        break;
      case Mesh::TETRAHEDRA:
        GeometricFactorsTet3D();
        break;
      case Mesh::HEXAHEDRA:
        GeometricFactorsHex3D();
        break;
    }
  }
  void GeometricFactorsTri2D();
  void GeometricFactorsTri3D();
  void GeometricFactorsQuad2D();
  void GeometricFactorsQuad3D();
  void GeometricFactorsTet3D();
  void GeometricFactorsHex3D();

  void SurfaceGeometricFactors() {
    switch (elementType) {
      case Mesh::TRIANGLES:
        if(dim==2)
          SurfaceGeometricFactorsTri2D();
        else
          SurfaceGeometricFactorsTri3D();
        break;
      case Mesh::QUADRILATERALS:
        if(dim==2)
          SurfaceGeometricFactorsQuad2D();
        else
          SurfaceGeometricFactorsQuad3D();
        break;
      case Mesh::TETRAHEDRA:
        SurfaceGeometricFactorsTet3D();
        break;
      case Mesh::HEXAHEDRA:
        SurfaceGeometricFactorsHex3D();
        break;
    }
  }
  void SurfaceGeometricFactorsTri2D();
  void SurfaceGeometricFactorsTri3D();
  void SurfaceGeometricFactorsQuad2D();
  void SurfaceGeometricFactorsQuad3D();
  void SurfaceGeometricFactorsTet3D();
  void SurfaceGeometricFactorsHex3D();

  void CubatureSetupTri2D();
  void CubatureSetupQuad2D();
  void CubatureSetupTet3D();
  void CubatureSetupHex3D();

  void CubaturePhysicalNodesTri2D();
  void CubaturePhysicalNodesTri3D();
  void CubaturePhysicalNodesQuad2D();
  void CubaturePhysicalNodesQuad3D();
  void CubaturePhysicalNodesTet3D();
  void CubaturePhysicalNodesHex3D();

  void PlotInterpTri2D(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch);
  void PlotInterpQuad2D(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch);
  void PlotInterpTet3D(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch);
  void PlotInterpHex3D(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch);

  void MassMatrixKernelSetupTri2D(int Nfields);
  void MassMatrixKernelSetupQuad2D(int Nfields);
  void MassMatrixKernelSetupTet3D(int Nfields);
  void MassMatrixKernelSetupHex3D(int Nfields);

  dfloat ElementCharacteristicLengthTri2D(dlong e);
  dfloat ElementCharacteristicLengthQuad2D(dlong e);
  dfloat ElementCharacteristicLengthTet3D(dlong e);
  dfloat ElementCharacteristicLengthHex3D(dlong e);

  /***************************************************************************/
  // Basic codes for generating nodes, polynomials, matrices, etc.

 public:
  //1D
  static void Nodes1D(const int _N, memory<dfloat>& _r);
  static void EquispacedNodes1D(const int _N, memory<dfloat>& _r);
  static void OrthonormalBasis1D(const dfloat a, const int i, dfloat& P);
  static void GradOrthonormalBasis1D(const dfloat a, const int i, dfloat& Pr);
  static void Vandermonde1D(const int _N,
                            const memory<dfloat> _r,
                            memory<dfloat>& V);
  static void GradVandermonde1D(const int _N,
                                const memory<dfloat> _r,
                                memory<dfloat>& Vr);

  static void MassMatrix1D(const int _Np,
                           const memory<dfloat> V,
                           memory<dfloat>& _MM);
  static void Dmatrix1D(const int _N,
                        const memory<dfloat> _rIn,
                        const memory<dfloat> _rOut,
                        memory<dfloat>& _Dr);
  static void InterpolationMatrix1D(const int _N,
                                    const memory<dfloat> _rIn,
                                    const memory<dfloat> _rOut,
                                    memory<dfloat>& I);
  static void DegreeRaiseMatrix1D(const int Nc, const int Nf,
                                  memory<dfloat>& P);
  static void CubatureWeakDmatrix1D(const int _Nq, const int _cubNq,
                                    const memory<dfloat> _cubProject,
                                    const memory<dfloat> _cubD,
                                    memory<dfloat>& _cubPDT);

  //Jacobi polynomial evaluation
  static dfloat JacobiP(const dfloat a, const dfloat alpha,
                        const dfloat beta, const int _N);
  static dfloat GradJacobiP(const dfloat a, const dfloat alpha,
                            const dfloat beta, const int _N);

  //Gauss-Legendre-Lobatto quadrature nodes
  static void JacobiGLL(const int _N,
                        memory<dfloat>& _x);
  static void JacobiGLL(const int _N,
                        memory<dfloat>& _x,
                        memory<dfloat>& _w);

  //Nth order Gauss-Jacobi quadrature nodes and weights
  static void JacobiGQ(const dfloat alpha, const dfloat beta,
                       const int _N,
                       memory<dfloat>& _x,
                       memory<dfloat>& _w);

  //Tris
  static void NodesTri2D(const int _N,
                         memory<dfloat>& _r,
                         memory<dfloat>& _s);
  static void FaceNodesTri2D(const int _N,
                             const memory<dfloat> _r,
                             const memory<dfloat> _s,
                             memory<int>& _faceNodes);
  static void VertexNodesTri2D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               memory<int>& _vertexNodes);
  static void FaceNodeMatchingTri2D(const memory<dfloat> _r,
                                    const memory<dfloat> _s,
                                    const memory<int> _faceNodes,
                                    const memory<int> _faceVertices,
                                    memory<int>& R);
  static void EquispacedNodesTri2D(const int _N,
                                   memory<dfloat>& _r,
                                   memory<dfloat>& _s);
  static void EquispacedEToVTri2D(const int _N, memory<int>& _EToV);
  static void SEMFEMNodesTri2D(const int _N,
                               int& _Np,
                               memory<dfloat>& _r,
                               memory<dfloat>& _s);
  static void SEMFEMEToVTri2D(const int _N,
                              int& _NelFEM,
                              memory<int>& _EToV);
  static void OrthonormalBasisTri2D(const dfloat _r, const dfloat _s,
                                    const int i, const int j,
                                    dfloat& P);
  static void GradOrthonormalBasisTri2D(const dfloat _r, const dfloat _s,
                                        const int i, const int j,
                                        dfloat& Pr, dfloat& Ps);
  static void VandermondeTri2D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               memory<dfloat>& V);
  static void GradVandermondeTri2D(const int _N,
                                   const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   memory<dfloat>& Vr,
                                   memory<dfloat>& Vs);
  static void MassMatrixTri2D(const int _Np,
                              const memory<dfloat> V,
                              memory<dfloat>& _MM);
  static void invMassMatrixTri2D(const int _Np,
                                 const memory<dfloat> V,
                                 memory<dfloat>& _invMM);
  static void DmatrixTri2D(const int _N,
                           const memory<dfloat> _r,
                           const memory<dfloat> _s,
                           memory<dfloat>& _D);
  static void LIFTmatrixTri2D(const int _N,
                              const memory<int> _faceNodes,
                              const memory<dfloat> _r,
                              const memory<dfloat> _s,
                              memory<dfloat>& _LIFT);
  static void SurfaceMassMatrixTri2D(const int _N,
                                     const memory<dfloat> _MM,
                                     const memory<dfloat> _LIFT,
                                     memory<dfloat>& _sM);
  static void SmatrixTri2D(const int _N,
                           const memory<dfloat> _Dr,
                           const memory<dfloat> _Ds,
                           const memory<dfloat> _MM,
                           memory<dfloat>& _S);
  static void InterpolationMatrixTri2D(const int _N,
                                       const memory<dfloat> rIn,
                                       const memory<dfloat> sIn,
                                       const memory<dfloat> rOut,
                                       const memory<dfloat> sOut,
                                       memory<dfloat>& I);
  static void DegreeRaiseMatrixTri2D(const int Nc, const int Nf,
                                     memory<dfloat>& P);
  static void CubatureNodesTri2D(const int cubTriN,
                                 int& _cubNp,
                                 memory<dfloat>& cubTrir,
                                 memory<dfloat>& cubTris,
                                 memory<dfloat>& cubTriw);
  static void CubaturePmatrixTri2D(const int _N,
                                   const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   const memory<dfloat> _cubr,
                                   const memory<dfloat> _cubs,
                                   memory<dfloat>& _cubProject);
  static void CubatureWeakDmatricesTri2D(const int _N,
                                         const memory<dfloat> _r,
                                         const memory<dfloat> _s,
                                         const memory<dfloat> _cubr,
                                         const memory<dfloat> _cubs,
                                         memory<dfloat>& _cubPDT);
  static void CubatureSurfaceMatricesTri2D(const int _N,
                                           const memory<dfloat> _r,
                                           const memory<dfloat> _s,
                                           const memory<int> _faceNodes,
                                           const memory<dfloat> _intr,
                                           const memory<dfloat> _intw,
                                           memory<dfloat>& _intInterp,
                                           memory<dfloat>& _intLIFT);
  static void SEMFEMInterpMatrixTri2D(const int _N,
                                      const memory<dfloat> _r,
                                      const memory<dfloat> _s,
                                      const memory<dfloat> _rFEM,
                                      const memory<dfloat> _sFEM,
                                      memory<dfloat>& I);

  static void Warpfactor(const int _N,
                         const memory<dfloat> _r,
                         memory<dfloat> warp);
  static void WarpBlendTransformTri2D(const int _N,
                                      memory<dfloat> _r,
                                      memory<dfloat> _s,
                                      const dfloat alphaIn=-1);


  //Quads
  static void NodesQuad2D(const int _N,
                          memory<dfloat>& _r,
                          memory<dfloat>& _s);
  static void FaceNodesQuad2D(const int _N,
                              const memory<dfloat> _r,
                              const memory<dfloat> _s,
                              memory<int>& _faceNodes);
  static void VertexNodesQuad2D(const int _N,
                                const memory<dfloat> _r,
                                const memory<dfloat> _s,
                                memory<int>& _vertexNodes);
  static void FaceNodeMatchingQuad2D(const memory<dfloat> _r,
                                     const memory<dfloat> _s,
                                     const memory<int> _faceNodes,
                                     const memory<int> _faceVertices,
                                     memory<int>& R);
  static void EquispacedNodesQuad2D(const int _N,
                                    memory<dfloat>& _r,
                                    memory<dfloat>& _s);
  static void EquispacedEToVQuad2D(const int _N, memory<int>& _EToV);
  static void SEMFEMEToVQuad2D(const int _N, memory<int>& _EToV);
  static void OrthonormalBasisQuad2D(const dfloat a, const dfloat b,
                                     const int i, const int j,
                                     dfloat& P);
  static void GradOrthonormalBasisQuad2D(const dfloat a, const dfloat b,
                                         const int i, const int j,
                                         dfloat& Pr, dfloat& Ps);
  static void VandermondeQuad2D(const int _N,
                                const memory<dfloat> _r,
                                const memory<dfloat> _s,
                                memory<dfloat>& V);
  static void GradVandermondeQuad2D(const int _N,
                                    const memory<dfloat> _r,
                                    const memory<dfloat> _s,
                                    memory<dfloat>& Vr,
                                    memory<dfloat>& Vs);
  static void MassMatrixQuad2D(const int _Np,
                               const memory<dfloat> V,
                               memory<dfloat>& _MM);
  static void LumpedMassMatrixQuad2D(const int _N,
                                     const memory<dfloat> _gllw,
                                     memory<dfloat>& _MM);
  static void invLumpedMassMatrixQuad2D(const int _N,
                                        const memory<dfloat> _gllw,
                                        memory<dfloat>& _invMM);
  static void DmatrixQuad2D(const int _N,
                            const memory<dfloat> _r,
                            const memory<dfloat> _s,
                            memory<dfloat>& _D);
  static void InterpolationMatrixQuad2D(const int _N,
                                        const memory<dfloat> rIn,
                                        const memory<dfloat> sIn,
                                        const memory<dfloat> rOut,
                                        const memory<dfloat> sOut,
                                        memory<dfloat>& I);

  //Tets
  static void NodesTet3D(const int _N,
                         memory<dfloat>& _r,
                         memory<dfloat>& _s,
                         memory<dfloat>& _t);
  static void FaceNodesTet3D(const int _N,
                             const memory<dfloat> _r,
                             const memory<dfloat> _s,
                             const memory<dfloat> _t,
                             memory<int>& _faceNodes);
  static void VertexNodesTet3D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               const memory<dfloat> _t,
                               memory<int>& _vertexNodes);
  static void FaceNodeMatchingTet3D(const memory<dfloat> _r,
                                    const memory<dfloat> _s,
                                    const memory<dfloat> _t,
                                    const memory<int> _faceNodes,
                                    const memory<int> _faceVertices,
                                    memory<int>& R);
  static void EquispacedNodesTet3D(const int _N,
                                   memory<dfloat>& _r,
                                   memory<dfloat>& _s,
                                   memory<dfloat>& _t);
  static void EquispacedEToVTet3D(const int _N, memory<int>& _EToV);
  static void SEMFEMEToVTet3D(const int _N, memory<int>& _EToV);
  static void OrthonormalBasisTet3D(const dfloat _r, const dfloat _s, const dfloat _t,
                                    const int i, const int j, const int k,
                                    dfloat& P);
  static void GradOrthonormalBasisTet3D(const dfloat _r, const dfloat _s, const dfloat _t,
                                        const int i, const int j, const int k,
                                        dfloat& Pr, dfloat& Ps, dfloat& Pt);
  static void VandermondeTet3D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               const memory<dfloat> _t,
                               memory<dfloat>& V);
  static void GradVandermondeTet3D(const int _N,
                                   const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   const memory<dfloat> _t,
                                   memory<dfloat>& Vr,
                                   memory<dfloat>& Vs,
                                   memory<dfloat>& Vt);
  static void MassMatrixTet3D(const int _Np,
                              const memory<dfloat> V,
                              memory<dfloat>& _MM);
  static void invMassMatrixTet3D(const int _Np,
                                 const memory<dfloat> V,
                                 memory<dfloat>& _invMM);
  static void DmatrixTet3D(const int _N,
                           const memory<dfloat> _r,
                           const memory<dfloat> _s,
                           const memory<dfloat> _t,
                           memory<dfloat>& _D);
  static void LIFTmatrixTet3D(const int _N,
                              const memory<int> _faceNodes,
                              const memory<dfloat> _r,
                              const memory<dfloat> _s,
                              const memory<dfloat> _t,
                              memory<dfloat>& _LIFT);
  static void SurfaceMassMatrixTet3D(const int _N,
                                     const memory<dfloat> _MM,
                                     const memory<dfloat> _LIFT,
                                     memory<dfloat>& _sM);
  static void SmatrixTet3D(const int _N,
                           const memory<dfloat> _Dr,
                           const memory<dfloat> _Ds,
                           const memory<dfloat> _Dt,
                           const memory<dfloat> _MM,
                           memory<dfloat>& _S);
  static void InterpolationMatrixTet3D(const int _N,
                                       const memory<dfloat> rIn,
                                       const memory<dfloat> sIn,
                                       const memory<dfloat> tIn,
                                       const memory<dfloat> rOut,
                                       const memory<dfloat> sOut,
                                       const memory<dfloat> tOut,
                                       memory<dfloat>& I);
  static void DegreeRaiseMatrixTet3D(const int Nc, const int Nf,
                                     memory<dfloat>& P);
  static void CubatureNodesTet3D(const int cubTetN,
                                 int& _cubNp,
                                 memory<dfloat>& _cubr,
                                 memory<dfloat>& _cubs,
                                 memory<dfloat>& _cubt,
                                 memory<dfloat>& _cubw);
  static void CubaturePmatrixTet3D(const int _N,
                                   const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   const memory<dfloat> _t,
                                   const memory<dfloat> _cubr,
                                   const memory<dfloat> _cubs,
                                   const memory<dfloat> _cubt,
                                   memory<dfloat>& _cubProject);
  static void CubatureWeakDmatricesTet3D(const int _N,
                                         const memory<dfloat> _r,
                                         const memory<dfloat> _s,
                                         const memory<dfloat> _t,
                                         const memory<dfloat> _cubr,
                                         const memory<dfloat> _cubs,
                                         const memory<dfloat> _cubt,
                                         memory<dfloat>& _cubPDT);
  static void CubatureSurfaceMatricesTet3D(const int _N,
                                           const memory<dfloat> _r,
                                           const memory<dfloat> _s,
                                           const memory<dfloat> _t,
                                           const memory<int> _faceNodes,
                                           const memory<dfloat> _intr,
                                           const memory<dfloat> _ints,
                                           const memory<dfloat> _intw,
                                           memory<dfloat>& _intInterp,
                                           memory<dfloat>& _intLIFT);
  static void SEMFEMInterpMatrixTet3D(const int _N,
                                      const memory<dfloat> _r,
                                      const memory<dfloat> _s,
                                      const memory<dfloat> _t,
                                      const memory<dfloat> _rFEM,
                                      const memory<dfloat> _sFEM,
                                      const memory<dfloat> _tFEM,
                                      memory<dfloat>& I);
  static void WarpShiftFace3D(const int _N, const dfloat alpha,
                              const memory<dfloat> L1,
                              const memory<dfloat> L2,
                              const memory<dfloat> L3,
                              memory<dfloat> w1,
                              memory<dfloat> w2);
  static void WarpBlendTransformTet3D(const int _N,
                                      memory<dfloat> _r,
                                      memory<dfloat> _s,
                                      memory<dfloat> _t,
                                      const dfloat alphaIn=-1);


  //Hexs
  static void NodesHex3D(const int _N,
                         memory<dfloat>& _r,
                         memory<dfloat>& _s,
                         memory<dfloat>& _t);
  static void FaceNodesHex3D(const int _N,
                             const memory<dfloat> _r,
                             const memory<dfloat> _s,
                             const memory<dfloat> _t,
                             memory<int>& _faceNodes);
  static void VertexNodesHex3D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               const memory<dfloat> _t,
                               memory<int>& _vertexNodes);
  static void FaceNodeMatchingHex3D(const memory<dfloat> _r,
                                    const memory<dfloat> _s,
                                    const memory<dfloat> _t,
                                    const memory<int> _faceNodes,
                                    const memory<int> _faceVertices,
                                    memory<int>& R);
  static void EquispacedNodesHex3D(const int _N,
                                   memory<dfloat>& _r,
                                   memory<dfloat>& _s,
                                   memory<dfloat>& _t);
  static void EquispacedEToVHex3D(const int _N, memory<int>& _EToV);
  static void SEMFEMEToVHex3D(const int _N, memory<int>& _EToV);
  static void OrthonormalBasisHex3D(const dfloat a, const dfloat b, const dfloat c,
                                    const int i, const int j, const int k,
                                    dfloat& P);
  static void GradOrthonormalBasisHex3D(const dfloat a, const dfloat b, const dfloat c,
                                        const int i, const int j, const int k,
                                        dfloat& Pr, dfloat& Ps, dfloat& Pt);
  static void VandermondeHex3D(const int _N,
                               const memory<dfloat> _r,
                               const memory<dfloat> _s,
                               const memory<dfloat> _t,
                               memory<dfloat>& V);
  static void GradVandermondeHex3D(const int _N,
                                   const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   const memory<dfloat> _t,
                                   memory<dfloat>& Vr,
                                   memory<dfloat>& Vs,
                                   memory<dfloat>& Vt);
  static void MassMatrixHex3D(const int _Np,
                              const memory<dfloat> V,
                              memory<dfloat>& _MM);
  static void LumpedMassMatrixHex3D(const int _N,
                                    const memory<dfloat> _gllw,
                                    memory<dfloat>& _MM);
  static void invLumpedMassMatrixHex3D(const int _N,
                                       const memory<dfloat> _gllw,
                                       memory<dfloat>& _invMM);
  static void DmatrixHex3D(const int _N,
                           const memory<dfloat> _r,
                           const memory<dfloat> _s,
                           const memory<dfloat> _t,
                           memory<dfloat>& _D);
  static void InterpolationMatrixHex3D(const int _N,
                                       const memory<dfloat> rIn,
                                       const memory<dfloat> sIn,
                                       const memory<dfloat> tIn,
                                       const memory<dfloat> rOut,
                                       const memory<dfloat> sOut,
                                       const memory<dfloat> tOut,
                                       memory<dfloat>& I);
};

} //namespace libp

#endif

