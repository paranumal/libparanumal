function computepoints(conn, mesh, r1d)
  N = length(r1d)-1

  r = repeat(r1d, outer=(1, N+1, N+1))
  s = repeat(r1d', outer=(N+1, 1, N+1))
  t = reshape(repeat(r1d, inner=((N+1)^2, 1, 1)), N+1, N+1, N+1)

  x = similar(r1d, N+1, N+1, N+1, mesh.Klocal)
  y = similar(r1d, N+1, N+1, N+1, mesh.Klocal)
  z = similar(r1d, N+1, N+1, N+1, mesh.Klocal)

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

  x, y, z
end
