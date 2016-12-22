
function writeNodeDataQuad2D(N)

r1d = JacobiGL(0,0,N);
V1d = Vandermonde1D(N, r1d);
D1d = Dmatrix1D(N, r1d, V1d);
M1d = inv(V1d')/V1d;
w1d = sum(M1d);
cnt = 1;
for i=1:N+1
  for j=1:N+1
    r(cnt) = r1d(j); %% r runs fastest
    s(cnt) = r1d(i);
    cnt = cnt+1;
  end
end
Np = (N+1)*(N+1);
r = reshape(r, Np,1);
s = reshape(s, Np,1);

Nfp = N+1;
Nfaces = 4;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)'; 
faceNodes2   = find( abs(r-1) < NODETOL)';
faceNodes3   = find( abs(s-1) < NODETOL)';
faceNodes4   = find( abs(r+1) < NODETOL)';
faceNodes  = [faceNodes1;faceNodes2;faceNodes3;faceNodes4]';

V = VandermondeQuad2D(N, r, s);
[Dr,Ds] = DmatricesQuad2D(N, r, s, V);
LIFT = LiftQuad2D(N, faceNodes, r, s);

fname = sprintf('quadrilateralN%02d.dat', N);

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

fprintf(fid, '%% faceNodes\n');
for f=1:Nfaces
  for m=1:Nfp
    fprintf(fid, '%d ', faceNodes(m,f)-1); %% adjust for 0-indexing
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
[plotR,plotS] = meshgrid(linspace(-1,1,N+4));
plotR = plotR(:); plotS = plotS(:);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunay(plotR,plotS)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = VandermondeQuad2D(N, plotR,plotS)/V; 

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


fclose(fid);
