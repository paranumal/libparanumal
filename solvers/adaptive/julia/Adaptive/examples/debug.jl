using MPI
using Pxest
using Pxest.p8est
using SparseArrays
using Adaptive
using LinearAlgebra

p8est_face_edges = [4 6  8 10
                    5 7  9 11
                    0 2  8  9
                    1 3 10 11
                    0 1  4  5
                    2 3  6  7]
p8est_corner_faces = [0 2 4
                      1 2 4
                      0 3 4
                      1 3 4
                      0 2 5
                      1 2 5
                      0 3 5
                      1 3 5]
p8est_corner_edges = [0 4  8
                      0 5  9
                      1 4 10
                      1 5 11
                      2 6  8
                      2 7  9
                      3 6 10
                      3 7 11]
p8est_corner_face_corners = [ 0 -1  0 -1  0 -1
                             -1  0  1 -1  1 -1
                              1 -1 -1  0  2 -1
                             -1  1 -1  1  3 -1
                              2 -1  2 -1 -1  0
                             -1  2  3 -1 -1  1
                              3 -1 -1  2 -1  2
                             -1  3 -1  3 -1  3]
p8est_edge_corners = [0  1
                      2  3
                      4  5
                      6  7
                      0  2
                      1  3
                      4  6
                      5  7
                      0  4
                      1  5
                      2  6
                      3  7]

function get_element_interpolation_operator(T, N, fc)
  r, w = lglpoints(T, N)
  Irb = interpolationmatrix(r, (r.-1)./2)
  Irt = interpolationmatrix(r, (r.+1)./2)

  #   hanging_face = -1 if the face is not hanging,
  #                = the corner of the full face that it touches:
  #                  e.g. if face = i and hanging_face[i] =
  #                  j, then the interpolation operator corresponding
  #                  to corner j should be used for that face.
  #  hanging_edge = -1 if the edge is not hanging,
  #               =  0 if the edge is the first half of a full edge,
  #                    but neither of the two faces touching the
  #                    edge is hanging,
  #               =  1 if the edge is the second half of a full edge,
  #                    but neither of the two faces touching the
  #                    edge is hanging,
  #               =  2 if the edge is the first half of a full edge
  #                    and is on the boundary of a full face,
  #               =  3 if the edge is the second half of a full edge
  #                    and is on the boundary of a full face,
  #               =  4 if the edge is in the middle of a full face.
  #                    See the diagram below for clarification.
  #
  # o...............o  o...............o  +---2---+.......o  o.......+---3---+
  # :               :  :               :  |       |       :  :       |       |
  # :               :  :               :  3   2   4       :  :       4   3   3
  # :               :  :               :  |       |       :  :       |       |
  # +---4---+       :  :       +---4---+  +---4---+       :  :       +---4---+
  # |       |       :  :       |       |  :               :  :               :
  # 2   0   4       :  :       4   1   2  :               :  :               :
  # |       |       :  :       |       |  :               :  :               :
  # +---2---+.......o  o.......+---3---+  o...............o  o...............o
  #
  #                    o                  +-------+
  #                    :                  |\       \
  #                    :                  1 \       \
  #                    :                  |  +-------+
  #                    +-------+          +  |       |
  #                    |\       \         :\ |       |
  #                    0 \       \        : \|       |
  #                    |  +-------+       :  +-------+
  #                    +  |       |       o
  #                     \ |       |
  #                      \|       |
  #                       +-------+
  #

  hf = -ones(Int, 6)
  he = -ones(Int, 12)
  c = fc & 0x0007
  work = fc >> 3
  cwork = c
  for i = 0:2
    if (work & 0x0001) != 0
      f = p8est_corner_faces[c+1,i+1]
      hf[f+1] = p8est_corner_face_corners[c+1,f+1]
      for j = 0:3
        e = p8est_face_edges[f+1,j+1]
        he[e+1] = 4
      end
    end
    work >>= 1;
  end
  for i = 0:2
    if (work & 0x0001) != 0
      e = p8est_corner_edges[c+1,i+1]
      he[e+1] = (he[e+1] == -1) ? 0 : 2
      he[e+1] += (cwork & 0x0001)
    end
    cwork >>= 1
    work >>= 1
  end

  II = zeros(T, (N+1)^3, (N+1)^3)

  l = LinearIndices((1:N+1, 1:N+1, 1:N+1))

  fmask = [l[1, :, :], l[N+1, :, :],
           l[:, 1, :], l[:, N+1, :],
           l[:, :, 1], l[:, :, N+1]]

  gmask = [l[  :,   1,   1],
           l[  :, N+1,   1],
           l[  :,   1, N+1],
           l[  :, N+1, N+1],
           l[  1,   :,   1],
           l[N+1,   :,   1],
           l[  1,   :, N+1],
           l[N+1,   :, N+1],
           l[  1,   1,   :],
           l[N+1,   1,   :],
           l[  1, N+1,   :],
           l[N+1, N+1,   :]]

  # Build up the interpolation matrix by applying the interpolations
  # to the identity matrix
  for i = 1:(N+1)^3
    q = zeros(T, N+1, N+1, N+1)
    q[i] = one(T)

    Iq = copy(q)
    for f = 0:5
      c = hf[f+1]

      if c > -1
        Ir1 = ((c & 1) == 0) ? Irb : Irt
        Ir2 = ((c & 2) == 0) ? Irb : Irt

        fm = reshape(fmask[f+1], :)
        fq = @view q[fm]

        Iq[fm] .= kron(Ir1, Ir2) * fq
      end
    end

    for g = 0:11
      c = he[g+1]

      if c > -1
        Ir1 = ((c & 1) == 0) ? Irb : Irt

        gm = reshape(gmask[g+1], :)
        gq = @view q[gm]

        Iq[gm] .= Ir1 * gq
      end
    end

    II[:,i] .= reshape(Iq,:)
  end

  II
end


function dump_vtk(pxest;
                  mpicomm = MPI.COMM_WORLD,
                  mpirank = MPI.Comm_rank(mpicomm),
                  vtk_dir = "vtk_files",
                  vtk_base = "mesh_p8est",
                 )
  # dump VTK
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  p8est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      force=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      force=true)
  end
end

!MPI.Initialized() && MPI.Init()

let
  N = 3 # polynomial degree
  refine_level = 1

  conn = p8est.Connectivity(1,1,2)
  pxest = p8est.PXEST(conn; min_lvl=0)

  p8est.refine!(pxest) do which_tree, quadrant
    qid = ccall(:p8est_quadrant_child_id, Cint,
                (Ref{Pxest.p8est.pxest_quadrant_t},),
                quadrant)
    add = (qid == 0 || qid == 3 || qid == 5 || qid == 6) ? 1 : 0

    refine =  quadrant.level < refine_level + add

    refine ? Cint(1) : Cint(0)
  end
  p8est.balance!(pxest)
  p8est.partition!(pxest)
  p8est.ghost!(pxest)
  p8est.lnodes!(pxest, N)
  mesh = p8est.Mesh(pxest)

  # WARNING assume single MPI Rank
  S = sparse(1:length(mesh.DToC), mesh.DToC[:], ones(Int, length(mesh.DToC)))
  G = S'

  r, w = lglpoints(Float64, N)
  D = spectralderivative(r)

  for fc in mesh.EToFC
    Ie = get_element_interpolation_operator(Float64, N, fc)

    # @show Ie
  end


  dump_vtk(pxest)
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  GC.gc()
  MPI.Finalize()
end
