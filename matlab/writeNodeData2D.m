
function writeNodeData2D(N)

[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);

Np = length(r);
Nfp = N+1;
Nfaces = 3;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)'; 
faceNodes2   = find( abs(r+s) < NODETOL)';
faceNodes3   = find( abs(r+1) < NODETOL)';
FaceNodes  = [faceNodes1;faceNodes2;faceNodes3]';

V = Vandermonde2D(N, r, s);
[Dr,Ds] = Dmatrices2D(N, r, s, V);
LIFT = Lift2D(N, FaceNodes, r, s);

fname = sprintf('triangleN%02d.dat', N);

fid = fopen(fname, 'w');

fprintf(fid, '%% degree N\n');
fprintf(fid, '%d\n', N);
fprintf(fid, '%% number of nodes\n');
fprintf(fid, '%d\n', Np);
fprintf(fid, '%% node coordinates\n');
for n=1:Np
  fprintf(fid, '%17.15E %17.15E\n', r(n), s(n));
end

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

fprintf(fid, '%% FaceNodes\n');
for f=1:Nfaces
  for m=1:Nfp
    fprintf(fid, '%d ', FaceNodes(m,f)-1); %% adjust for 0-indexing
  end
  fprintf(fid, '\n');
end

fprintf(fid, '%% LIFT matrix\n');
for n=1:Np
  for m=1:Nfp*Nfaces
    fprintf(fid, '%17.15E ', LIFT(n,m));
  end
  fprintf(fid, '\n');
end

%% compute equispaced nodes on equilateral triangle
[plotR,plotS] = EquiNodes2D(N+4);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunay(plotR,plotS)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%% transform to bi-unit triangle
[plotR,plotS] = xytors(plotR,plotS);

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = Vandermonde2D(N, plotR,plotS)/V; 

%% output plot nodes
fprintf(fid, '%% number of plot nodes\n');
fprintf(fid, '%d\n', plotNp);
fprintf(fid, '%% plot node coordinates\n');
for n=1:plotNp
  fprintf(fid, '%17.15E %17.15E\n', plotR(n), plotS(n));
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
  fprintf(fid, '%d %d %d\n' ,...
 	plotEToV(n,1),plotEToV(n,2),plotEToV(n,3));
end

%% volume cubature
[cubr,cubs,cubw] = Cubature2D(3*N);
cInterp = Vandermonde2D(N, cubr, cubs)/V;
Np
Ncub = length(cubr)

fprintf(fid, '%% number of volume cubature nodes\n');
fprintf(fid, '%d\n', length(cubr));

fprintf(fid, '%% cubature node coordinates\n');
for n=1:Ncub
  fprintf(fid, '%17.15E %17.15E\n', cubr(n), cubs(n));
end

fprintf(fid, '%% cubature interpolation matrix\n');
for n=1:Ncub
  for m=1:Np
    fprintf(fid, '%17.15E ', cInterp(n,m));
  end
  fprintf(fid, '\n');
end

cV = Vandermonde2D(N, cubr, cubs);
cV'*diag(cubw)*cV;

[cVr,cVs] = GradVandermonde2D(N, cubr, cubs);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
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

fprintf(fid, '%% cubature projection matrix\n');
for n=1:Np
  for m=1:Ncub
    fprintf(fid, '%17.15E ', cubProject(n,m));
  end
  fprintf(fid, '\n');
end
cubProject
testIdentity = cubProject*cInterp;


%% surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
%z = JacobiGL(0,0,N);
%zV = Vandermonde1D(N,z);
%w = sum(inv(zV*transpose(zV)));

Nfi = length(z);

ir = [z,-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, -z];
iw = [w,w,w];

sV = Vandermonde2D(N, ir(:), is(:));
sInterp = sV/V;

iInterp = [sInterp(1:Nfi,FaceNodes(:,1));sInterp(Nfi+1:2*Nfi,FaceNodes(:,2));sInterp(2*Nfi+1:3*Nfi,FaceNodes(:,3))];

fprintf(fid, '%% number of surface integration nodes per face\n');
fprintf(fid, '%d\n', length(z));
fprintf(fid, '%% surface integration interpolation matrix\n');
for n=1:Nfi*Nfaces
  for m=1:Nfp
    fprintf(fid, '%17.15E ', iInterp(n,m));
  end
  fprintf(fid, '\n');
end

bInterp = [];
bInterp(1:Nfi,1:Nfp) = iInterp(1:Nfi,:);
bInterp(Nfi+1:2*Nfi,Nfp+1:2*Nfp) = iInterp(Nfi+1:2*Nfi,:);
bInterp(2*Nfi+1:3*Nfi,2*Nfp+1:3*Nfp) = iInterp(2*Nfi+1:3*Nfi,:);

%% integration node lift matrix
iLIFT = V*V'*sInterp'*diag(iw(:));
size(iLIFT)
size(iInterp)
max(max(abs(iLIFT*bInterp-LIFT)))

fprintf(fid, '%% surface integration lift matrix\n');
for n=1:Np
  for m=1:Nfi*Nfaces
    fprintf(fid, '%17.15E ', iLIFT(n,m));
  end
  fprintf(fid, '\n');
end

cubDrT*ones(Ncub,1) 
nr = [zeros(Nfi,1);ones(Nfi,1);-ones(Nfi,1)];
ns = [-ones(Nfi,1);ones(Nfi,1);zeros(Nfi,1)];
sJ = [ones(Nfi,1);ones(Nfi,1);ones(Nfi,1)];
cubDrT*ones(Ncub,1) - iLIFT*(nr.*sJ)
cubDsT*ones(Ncub,1) - iLIFT*(ns.*sJ)






fclose(fid)
