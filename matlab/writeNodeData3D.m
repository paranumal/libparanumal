
function writeNodeData3D(N)

[r,s,t] = Nodes3D(N);
[r,s,t] = xyztorst(r,s,t);

Np = length(r);
Nfp = (N+1)*(N+2)/2;
Nfaces = 4;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(t+1) < NODETOL)'; 
faceNodes2   = find( abs(s+1) < NODETOL)';
faceNodes3   = find( abs(r+s+t+1) < NODETOL)';
faceNodes4   = find( abs(r+1) < NODETOL)';
faceNodes  = [faceNodes1;faceNodes2;faceNodes3;faceNodes4]';

V = Vandermonde3D(N, r, s, t);
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
LIFT = Lift3D(N, faceNodes, r, s, t);

fname = sprintf('tetN%02d.dat', N);

fid = fopen(fname, 'w');

fprintf(fid, '%% degree N\n');
fprintf(fid, '%d\n', N);
fprintf(fid, '%% number of nodes\n');
fprintf(fid, '%d\n', Np);
fprintf(fid, '%% node coordinates\n');
for n=1:Np
  fprintf(fid, '%17.15E %17.15E %17.15E\n', r(n), s(n), t(n));
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

fprintf(fid, '%% t collocation differentation matrix\n');
for n=1:Np
  for m=1:Np
    fprintf(fid, '%17.15E ', Dt(n,m));
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

%% compute equispaced nodes on equilateral triangle
[plotR,plotS,plotT] = EquiNodes3D(N+4);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
%plotEToV = delaunay3(plotR,plotS,plotT)-1; 
plotEToV = delaunayFixVolume(plotR,plotS,plotT)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = Vandermonde3D(N, plotR,plotS,plotT)/V; 

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

%% output cubature arrays

if N < 7
    [cubr cubs cubt cubw] = tet_cubature(2*N+1);
    Vq = Vandermonde3D(N,cubr,cubs,cubt)/V;
    Pq = V*V'*Vq'*diag(cubw);
    
    Ncub = length(cubw);
    fprintf(fid, '%% number of volume cubature nodes\n');
    fprintf(fid, '%d\n', Ncub);
    
    fprintf(fid, '%% cubature node coordinates\n');
    for n=1:Ncub
        fprintf(fid, '%17.15E %17.15E %17.15E\n', cubr(n), cubs(n), cubt(n));
    end
    
    fprintf(fid,'%% cubature interpolation matrix\n');
    for i=1:Ncub
        for j = 1:Np
            fprintf(fid, '%17.15E ', Vq(i,j));
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid,'%% cubature projection matrix\n');
    for i=1:Np
        for j = 1:Ncub
            fprintf(fid, '%17.15E ', Pq(i,j));
        end
        fprintf(fid, '\n');
    end
end
%%

fclose(fid);
