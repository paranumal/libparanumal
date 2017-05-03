
function writeNodeData3D(inN)

Globals3D;
N = inN;
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

MM = inv(transpose(V))/V;
fprintf(fid, '%% reference mass matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', MM(n,m));
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
    
    cV = Vandermonde3D(N, cubr, cubs, cubt);
    [cVr,cVs,cVt] = GradVandermonde3D(N, cubr, cubs, cubt);
    cubDrT = V*transpose(cVr)*diag(cubw);
    cubDsT = V*transpose(cVs)*diag(cubw);
    cubDtT = V*transpose(cVt)*diag(cubw);

    Ncub = length(cubw);
    fprintf(fid, '%% number of volume cubature nodes\n');
    fprintf(fid, '%d\n', Ncub);
    
    fprintf(fid, '%% cubature node coordinates\n');
    for n=1:Ncub
        fprintf(fid, '%17.15E %17.15E %17.15E %17.15E\n', cubr(n), cubs(n), cubt(n), cubw(n));
    end
    
    fprintf(fid,'%% cubature interpolation matrix\n');
    for i=1:Ncub
        for j = 1:Np
            fprintf(fid, '%17.15E ', Vq(i,j));
        end
        fprintf(fid, '\n');
    end

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

    
    fprintf(fid,'%% cubature projection matrix\n');
    for i=1:Np
        for j = 1:Ncub
            fprintf(fid, '%17.15E ', Pq(i,j));
        end
        fprintf(fid, '\n');
    end


end
%%

%% elliptic patch problem
K = 5;

VX = [-1, 1,      0,          0,           0,         5/3,        -5/3,         0];
VY = [ 0, 0,sqrt(3),  1/sqrt(3),-7*sqrt(3)/9, 8*sqrt(3)/9, 8*sqrt(3)/9, 1/sqrt(3)];
VZ = [ 0, 0,      0,2*sqrt(6)/3, 4*sqrt(6)/9, 4*sqrt(6)/9, 4*sqrt(6)/9,-2*sqrt(6)/3];

EToV = [1,2,3,4;
        1,2,4,5;
        2,3,4,6;
        3,1,4,7;
        1,3,2,8];

BCType = [0,0,0,0;
          0,0,0,0;
          0,0,0,0;
          0,0,0,0;
          0,0,0,0];

StartUp3D;
  
% build weak Poisson operator matrices
[A, M] = PoissonIPDG3D();

%% hack since we know W&B face 1 nodes are first
vmapP = reshape(vmapP, Nfp*Nfaces, K);
idsP = vmapP(:,1);
subind = [(1:Np)';idsP];

subA = full(A(subind,subind));
subM = full(M(subind,subind));

spy(abs(subA)>1e-10);
condSubA = cond(subA);

[B,d] = eig(subA, subM);

%% A = S + lambda*M 
%%   = M*(M\S + lambda*I) 
%%   ~ J*Mhat*(Mhat\Shat/hscale2 + lambda*I) 
%%   ~ J*Mhat*Bhat*(d/scale + lambda*I)/Bhat
%% inv(A) ~ Bhat*inv(J*(d/scale+lambda*I))*inv(Mhat*Bhat)

forwardMatrix = inv(subM*B);
diagOp = diag(d);
backwardMatrix = B;

NpP = size(subA,1);
fprintf(fid, '%% stencil size for IPDG OAS NpP\n');
fprintf(fid, '%d\n', NpP);
fprintf(fid, '%% forward matrix [ change of basis for IPDG OAS precon ]\n');
for n=1:NpP
  for m=1:NpP
  fprintf(fid, '%17.15E ', forwardMatrix(n,m));
  end
  fprintf(fid, '\n');
end
fprintf(fid, '%% diagOp \n');
for n=1:NpP
    fprintf(fid, '%17.15E ', diagOp(n));
  fprintf(fid, '\n');
end
fprintf(fid, '%% backward matrix [ reverse change of basis for IPDG OAS precon ]\n');
for n=1:NpP
  for m=1:NpP
    fprintf(fid, '%17.15E ', backwardMatrix(n,m));
  end
  fprintf(fid, '\n');
end

fclose(fid);
