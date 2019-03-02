using MPI
using Pxest
using Pxest.p8est
using SparseArrays
using Adaptive
using LinearAlgebra

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

  Snc = spzeros(size(S,1),size(S,1))

  for i in 1:mesh.Klocal
    Ie = nonconinterpolation(Float64, N, mesh.EToFC[i])

    idx = (i-1)*(N+1)^3 .+ (1:(N+1)^3)
    Snc[idx, idx] .= Ie
  end

  Gnc = Snc'

  GS = G * Gnc * Snc * S

  r1d, w1d = lglpoints(Float64, N)

  r = repeat(r1d, outer=(1, N+1, N+1))
  s = repeat(r1d', outer=(N+1, 1, N+1))
  t = reshape(repeat(r1d, inner=((N+1)^2, 1, 1)), N+1, N+1, N+1)

  x = similar(Snc, N+1, N+1, N+1, mesh.Klocal)
  y = similar(Snc, N+1, N+1, N+1, mesh.Klocal)
  z = similar(Snc, N+1, N+1, N+1, mesh.Klocal)

  for e in 1:mesh.Klocal
    tree = mesh.EToT[e]

    v = conn.tree_to_vertex[:,tree+1] .+ 1

    cr = mesh.EToX[e] / Pxest.p8est.PXEST_ROOT_LEN
    cs = mesh.EToY[e] / Pxest.p8est.PXEST_ROOT_LEN
    ct = mesh.EToZ[e] / Pxest.p8est.PXEST_ROOT_LEN
    level = mesh.EToL[e]
    h = 1 / (1 << (level + 1))

    r = cr .+ h .* (r1d .+ 1)
    s = cs .+ h .* (r1d .+ 1)
    t = ct .+ h .* (r1d .+ 1)

    for k=1:N+1, j=1:N+1, i=1:N+1
      w0 = (1 - r[i]) * (1 - s[j]) * (1 - t[k])
      w1 = r[i] * (1 - s[j]) * (1 - t[k])
      w2 = (1 - r[i]) * s[j] * (1 - t[k])
      w3 = r[i] * s[j] * (1 - t[k])
      w4 = (1 - r[i]) * (1 - s[j]) * t[k]
      w5 = r[i] * (1 - s[j]) * t[k]
      w6 = (1 - r[i]) * s[j] * t[k]
      w7 = r[i] * s[j] * t[k]

      x[i,j,k,e] = w0 * conn.vertices[1,v[1]] + w1 * conn.vertices[1,v[2]] +
                   w2 * conn.vertices[1,v[3]] + w3 * conn.vertices[1,v[4]] +
                   w4 * conn.vertices[1,v[5]] + w5 * conn.vertices[1,v[6]] +
                   w6 * conn.vertices[1,v[7]] + w7 * conn.vertices[1,v[8]]

      y[i,j,k,e] = w0 * conn.vertices[2,v[1]] + w1 * conn.vertices[2,v[2]] +
                   w2 * conn.vertices[2,v[3]] + w3 * conn.vertices[2,v[4]] +
                   w4 * conn.vertices[2,v[5]] + w5 * conn.vertices[2,v[6]] +
                   w6 * conn.vertices[2,v[7]] + w7 * conn.vertices[2,v[8]]

      z[i,j,k,e] = w0 * conn.vertices[3,v[1]] + w1 * conn.vertices[3,v[2]] +
                   w2 * conn.vertices[3,v[3]] + w3 * conn.vertices[3,v[4]] +
                   w4 * conn.vertices[3,v[5]] + w5 * conn.vertices[3,v[6]] +
                   w6 * conn.vertices[3,v[7]] + w7 * conn.vertices[3,v[8]]
    end
  end

  maximum(x)

  D = spectralderivative(r)

  sJ = similar(x, (N+1)^2, 6, mesh.Klocal)
  nx = similar(x, (N+1)^2, 6, mesh.Klocal)
  ny = similar(x, (N+1)^2, 6, mesh.Klocal)
  nz = similar(x, (N+1)^2, 6, mesh.Klocal)

  rx, ry, rz = similar(x), similar(x), similar(x)
  sx, sy, sz = similar(x), similar(x), similar(x)
  tx, ty, tz = similar(x), similar(x), similar(x)

  J = similar(x)
  W = reshape(repeat(kron(1, w1d, w1d, w1d), outer=(mesh.Klocal,)),
              (N+1, N+1, N+1, mesh.Klocal))

  computemetric!(x, y, z, J, rx, sx, tx, ry, sy, ty, rz, sz, tz, sJ,
                 nx, ny, nz, D)

  writemesh("mesh", x, y, z)

  dump_vtk(pxest)
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  GC.gc()
  MPI.Finalize()
end
