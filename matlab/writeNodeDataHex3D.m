
function writeNodeDataHex3D(N)

Np = (N+1)*(N+1)*(N+1);

r1d = JacobiGL(0,0,N);
V1d = Vandermonde1D(N, r1d);
D1d = Dmatrix1D(N, r1d, V1d);
M1d = inv(V1d')/V1d;
w1d = sum(M1d);
cnt = 1;
r = zeros(Np,1);
s = zeros(Np,1);
t = zeros(Np,1);
for k=1:N+1
  for j=1:N+1
    for i=1:N+1
      r(cnt) = r1d(i); %% r runs fastest
      s(cnt) = r1d(j);
      t(cnt) = r1d(k);
      cnt = cnt+1;
    end
  end
end
Np = (N+1)*(N+1)*(N+1);
r = reshape(r, Np,1);
s = reshape(s, Np,1);
t = reshape(t, Np,1);

Nfp = (N+1)*(N+1);
Nfaces = 6;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes = zeros(Nfp, Nfaces);
faceNodes(:,1)   = find( abs(t+1) < NODETOL); 
faceNodes(:,2)   = find( abs(s+1) < NODETOL);
faceNodes(:,3)   = find( abs(r-1) < NODETOL);
faceNodes(:,4)   = find( abs(s-1) < NODETOL);
faceNodes(:,5)   = find( abs(r+1) < NODETOL);
faceNodes(:,6)   = find( abs(t-1) < NODETOL);

V = VandermondeHex3D(N, r, s, t);
[Dr,Ds,Dt] = DmatricesHex3D(N, r, s, t, V);
LIFT = LiftHex3D(N, faceNodes, r, s, t);

fname = sprintf('hexN%02d.dat', N);

fid = fopen(fname, 'w');

fprintf(fid, '%% degree N\n');
fprintf(fid, '%d\n', N);
fprintf(fid, '%% number of nodes\n');
fprintf(fid, '%d\n', Np);
fprintf(fid, '%% node coordinates\n');
for n=1:Np
  fprintf(fid, '%17.15E %17.15E %17.15E\n', r(n), s(n), t(n));
end

if(0)

fprintf(fid, '%% r collocation differentation matrix\n');
for n=1:Np
  for m=1:Np
    fprintf(fid, '%17.15E ', Dr(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% s collocation differentation matrix\n');
for n=1:Np
  for m=1:Np
    fprintf(fid, '%17.15E ', Ds(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% t collocation differentation matrix\n');
for n=1:Np
  for m=1:Np
    fprintf(fid, '%17.15E ', Dt(n,m));
  end
  fprintf(fid, '\n');
end
end

fprintf(fid, '%% faceNodes\n');
for f=1:Nfaces
  for m=1:Nfp
    fprintf(fid, '%d ', faceNodes(m,f)-1); %% adjust for 0-indexing
  end
  fprintf(fid, '\n');
end

if(0)
fprintf(fid, '%% LIFT matrix\n');
for n=1:Np
  for m=1:Nfp*Nfaces
    fprintf(fid, '%17.15E ', LIFT(n,m));
  end
  fprintf(fid, '\n');
end

end

fprintf(fid, '%% D (1D) matrix\n');
for n=1:N+1
  for m=1:N+1
    fprintf(fid, '%17.15E ', D1d(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% r (1D) gll nodes\n');
for n=1:N+1
fprintf(fid, '%17.15E \n', r1d(n));
end

fprintf(fid, '%% w (1D) gll node weights\n');
for n=1:N+1
fprintf(fid, '%17.15E \n', w1d(n));
end


%% compute equispaced nodes on equilateral triangle
[plotR,plotS,plotT] = meshgrid(linspace(-1,1,2)); %% hacked
plotR = plotR(:); plotS = plotS(:); plotT = plotT(:);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunay(plotR,plotS,plotT)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = VandermondeHex3D(N, plotR,plotS,plotT)/V; 

%% output plot nodes
fprintf(fid, '%% number of plot nodes\n');
fprintf(fid, '%d\n', plotNp);
fprintf(fid, '%% plot node coordinates\n');
for n=1:plotNp
  fprintf(fid, '%17.15E %17.15E %17.15E\n', plotR(n), plotS(n), plotT(n));
end

%% output plot interpolation matrix
fprintf(fid, '%% plot node interpolation matrix\n');
for n=1:plotNp
  for m=1:Np
    fprintf(fid, '%17.15E ', plotInterp(n,m));
  end
  fprintf(fid, '\n');
end

%% output plot triangulation
fprintf(fid, '%% number of plot elements\n');
fprintf(fid, '%d\n', plotNelements);

fprintf(fid, '%% number of vertices per plot elements\n');
fprintf(fid, '%d\n', size(plotEToV,2));

fprintf(fid, '%% triangulation of plot nodes\n');
for n=1:plotNelements
  fprintf(fid, '%d %d %d %d\n' ,...
 	plotEToV(n,1),plotEToV(n,2),plotEToV(n,3),plotEToV(n,4));
end

%%%% ---

%% 1D 
gllS = transpose(D1d)*diag(w1d)*D1d;

NqP = N+3;

%ids 
Nelements = 3;
cnt = 1;
for e=1:Nelements
for n=1:N+1
galnums(n,e) = cnt;
cnt = cnt+1;
end
cnt = cnt-1;
end

A = zeros(cnt,cnt);
for e=1:Nelements
for n=1:N+1
for m=1:N+1
i = galnums(n,e);
j = galnums(m,e);

A(i,j) = A(i,j) + gllS(n,m);
end
end
end

%% WARNING NEED N>1 (otherwise we need a boundary condition)

overlap = 1;
ids = N+1-overlap:2*N+1+overlap;
subA = A(ids,ids)

SP = zeros(NqP,NqP); %% one point overlap
SP(2:NqP-1,2:NqP-1) = gllS;
SP(1,1) = SP(1,1) + gllS(2,2);
SP(2,1) = SP(2,1) + gllS(1,2);
SP(2,2) = SP(2,2) + gllS(1,1);
SP(1,2) = SP(1,2) + gllS(2,1);

SP(NqP,NqP)   = SP(NqP,NqP) + gllS(2,2);
SP(NqP-1,NqP) = SP(NqP-1,NqP) + gllS(1,2);
SP(NqP-1,NqP-1) = SP(NqP-1,NqP-1) + gllS(1,1);
SP(NqP,NqP-1) = SP(NqP,NqP-1) + gllS(2,1);

gllwP = diag([w1d(2),2*w1d(1),w1d(2:end-1),2*w1d(1),w1d(2)]);

[vSP,dSP] = eig(gllwP\SP);

%% invSP = vSP*inv(dSP)*inv(gllwP*vSP) (i.e. the inverse we need)
%% define P = vSP, invP = inv(gllwP*vSP), 
P = vSP;
invP = inv(gllwP*vSP);
diagOp = diag(dSP); % need to divide to get inverse

fprintf(fid, '%% stencil size for OAS NqP\n');
fprintf(fid, '%d\n', NqP);
fprintf(fid, '%% invP [ change of basis for OAS precon ]\n');
for n=1:NqP
  for m=1:NqP
    fprintf(fid, '%17.15E ', invP(n,m));
  end
  fprintf(fid, '\n');
end
fprintf(fid, '%% diagOp [ weight so that inv(inv(W)*(trans(D)*W*D)) = P*inv(diagOp)*invP ]\n');
for n=1:NqP
    fprintf(fid, '%17.15E ', diagOp(n));
  fprintf(fid, '\n');
end
fprintf(fid, '%% P [ reverse change of basis for OAS precon ]\n');
for n=1:NqP
  for m=1:NqP
    fprintf(fid, '%17.15E ', P(n,m));
  end
  fprintf(fid, '\n');
end


%%ids 
Nelements = 10;
Nq = N+1;
ADG = zeros(Nq*Nelements,Nq*Nelements);
es = reshape(1:Nelements*Nq, Nq,Nelements);

tau = 2*(N+1)^2;
for e=2:Nelements-1
  n = es(:,e);
  nL = es(1,e);  
  nR = es(Nq,e);
  nP = es(:,e+1);
  nM = es(:,e-1);

  ADG(n,n)  = ADG(n,n)+gllS;

  ADG(n,nL)   = ADG(n,nL)   + 0.5*transpose(D1d(1,:));
  ADG(n,nL-1) = ADG(n,nL-1) - 0.5*transpose(D1d(1,:));

  ADG(n,nR)   = ADG(n,nR)   - 0.5*transpose(D1d(Nq,:));
  ADG(n,nR+1) = ADG(n,nR+1) + 0.5*transpose(D1d(Nq,:));

  ADG(nL,n)   = ADG(nL,n)   + 0.5*(D1d(1,:));
  ADG(nL-1,n) = ADG(nL-1,n) - 0.5*(D1d(1,:));

  ADG(nR,n)   = ADG(nR,n)   - 0.5*(D1d(Nq,:));
  ADG(nR+1,n) = ADG(nR+1,n) + 0.5*(D1d(Nq,:));

  ADG(nL,nL)  = ADG(nL,nL) + 0.5*tau;	    
  ADG(nL,nL-1) = ADG(nL,nL-1) - 0.5*tau;	    

  ADG(nR,nR) = ADG(nR,nR) + 0.5*tau;	    
  ADG(nR,nR+1) = ADG(nR,nR+1) - 0.5*tau;	    

  MDG(n,n) = diag(w1d);
end

ids = 4*Nq:5*Nq+1;
BDG = ADG(ids,ids)
MDG = MDG(ids,ids)

gllwP = diag([w1d(1),w1d,w1d(1)]);

[vSP,dSP] = eig(gllwP\BDG);

%% invSP = vSP*inv(dSP)*inv(gllwP*vSP) (i.e. the inverse we need)
%% define P = vSP, invP = inv(gllwP*vSP), 
P = vSP;
invP = inv(gllwP*vSP);
diagOp = diag(dSP); % need to divide to get inverse

fprintf(fid, '%% stencil size for DG OAS NqP\n');
fprintf(fid, '%d\n', NqP);
fprintf(fid, '%% forwardDG [ change of basis for DG OAS precon ]\n');
for n=1:NqP
  for m=1:NqP
    fprintf(fid, '%17.15E ', invP(n,m));
  end
  fprintf(fid, '\n');
end
fprintf(fid, '%% diagOpDG [ weight so that inv(inv(W)*(trans(D)*W*D)) = P*inv(diagOp)*invP ]\n');
for n=1:NqP
    fprintf(fid, '%17.15E ', diagOp(n));
  fprintf(fid, '\n');
end
fprintf(fid, '%% backwardDG [ reverse change of basis forr H0 OAS precon ]\n');
for n=1:NqP
  for m=1:NqP
    fprintf(fid, '%17.15E ', P(n,m));
  end
  fprintf(fid, '\n');
end


[gr,gw] = JacobiGQ(0,0,N+1);
gI = Vandermonde1D(N, gr)/Vandermonde1D(N, r1d);
gD = GradVandermonde1D(N, gr)/Vandermonde1D(N, r1d);
gNq = length(gr);

fprintf(fid, '%% gNq (gauss quadrature count)\n');
fprintf(fid, '%d\n', gNq);
fprintf(fid, '%% gjr,gjw Gauss Jacobi nodes and weights\n');
for n=1:gNq
  fprintf(fid, '%17.15E %17.15E \n', gr(n), gw(n));
end
fprintf(fid, '%% gjI [ 1D interpolation from GLL to Gauss nodes  ]\n');
for n=1:gNq
  for m=1:N+1
    fprintf(fid, '%17.15E ', gI(n,m));
  end
  fprintf(fid, '\n');
end
fprintf(fid, '%% gjD [ 1D differentiation from GLL to Gauss nodes  ]\n');
for n=1:gNq
  for m=1:N+1
    fprintf(fid, '%17.15E ', gD(n,m));
  end
  fprintf(fid, '\n');
end

fclose(fid);
if(0)


%% volume cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
Nz = length(z)
sk = 1;
cubr = zeros(Nz*Nz*Nz,1);
cubs = zeros(Nz*Nz*Nz,1);
cubt = zeros(Nz*Nz*Nz,1);
cubw = zeros(Nz*Nz*Nz,1);
for k=1:Nz
for j=1:Nz
for i=1:Nz
cubr(sk) = z(i);
cubs(sk) = z(j);
cubt(sk) = z(k);
cubw(sk) = w(i)*w(j)*w(k);
sk = sk+1;
end
end
end

cInterp = VandermondeHex3D(N, cubr, cubs, cubt)/V;
Ncub = length(cubr);

fprintf(fid, '%% number of volume cubature nodes\n');
fprintf(fid, '%d\n', length(cubr));
fprintf(fid, '%% cubature interpolation matrix\n');
for n=1:Ncub
  for m=1:Np
    fprintf(fid, '%17.15E ', cInterp(n,m));
  end
  fprintf(fid, '\n');
end

cV = VandermondeHex3D(N, cubr, cubs, cubt);
cV'*diag(cubw)*cV;

[cVr,cVs,cVt] = GradVandermondeHex3D(N, cubr, cubs, cubt);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
cubDtT = V*transpose(cVt)*diag(cubw);
cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)

fprintf(fid, '%% cubature r weak differentiation matrix\n');
for n=1:Np
  for m=1:Ncub
    fprintf(fid, '%17.15E ', cubDrT(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% cubature s weak differentiation matrix\n');
for n=1:Np
  for m=1:Ncub
    fprintf(fid, '%17.15E ', cubDsT(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% cubature t weak differentiation matrix\n');
for n=1:Np
  for m=1:Ncub
    fprintf(fid, '%17.15E ', cubDtT(n,m));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% cubature projection matrix\n');
for n=1:Np
  for m=1:Ncub
    fprintf(fid, '%17.15E ', cubProject(n,m));
  end
  fprintf(fid, '\n');
end
cubProject;
testIdentity = cubProject*cInterp;

fclose(fid);
end
