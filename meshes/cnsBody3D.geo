
// AGARD_R = 1, AGARD_O = 2, STANDALONE = 3
MODEL = 1; 

BBOX  = 1; // SPHERE 
// BBOX  = 2; // BOX 
// Note that Length = 8.5 x D for AGARD B
LoD = 8.5; 

SetFactory("OpenCASCADE");


If(MODEL == 1)
  v() = ShapeFromFile("agard_b_1.stp");
EndIf


If(MODEL == 2)
  v() = ShapeFromFile("agard_r.stp");
EndIf

If(MODEL == 3)
  v() = ShapeFromFile("bodyalone.stp");
EndIf


// Get the bounding box of the volume:
bbox() = BoundingBox Volume{v()};
xmin = bbox(0); ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3); ymax = bbox(4);zmax = bbox(5);
// Report the bounding box
For i In {0:#bbox()-1}
  Printf ("bounding box = %g", bbox(i));
EndFor

// Find the current center
xc = 0.5*(xmin + xmax);
yc = 0.5*(ymin + ymax);
zc = 0.5*(zmin + zmax);

scale_factor = LoD/(xmax - xmin); 

// Scale the volume so that (xmax - xmin) = 8.5; 
Dilate {{0, 0, 0}, {scale_factor, scale_factor, scale_factor}} { Volume{1};}

// Get the bounding box of the volume again:
bbox() = BoundingBox Volume{v()};
xmin = bbox(0);ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3);ymax = bbox(4);zmax = bbox(5);

// Get the center to [0 0 0]
xc = 0.5*(xmin + xmax); yc = 0.5*(ymin + ymax); zc = 0.5*(zmin + zmax);
Translate{-xc, -yc, -zc}{Volume{v()};}

bbox() = BoundingBox Volume{1};
xmin = bbox(0);ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3);ymax = bbox(4);zmax = bbox(5);


For i In {0:#bbox()-1}
  Printf ("bounding box after centering = %g", bbox(i));
EndFor

If(BBOX == 1)
  RI = 10*(xmax - xmin)/2; 
  Sphere(2) = {1, 0, 0, RI, -Pi/2, Pi/2, 2*Pi};
  BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
ElseIf(BBOX == 2)
  xbmin = -20; DX = 50;  
  ybmin = -15; DY = 30; 
  zbmin = -15; DZ = 30; 
  Box(2) = {xbmin, ybmin, zbmin, DX, DY, DZ};
  BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
EndIf
 

// Print out wall surface ids
vol()     = Volume In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
walls()   = Surface In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};
outSurf() = Surface In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
wpoints() = Point In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};

// Volume ID check
For i In {0:#vol()-1}
  Printf ("volume Cehck : id of # volume = %g / %g", vol(i),#vol() );
EndFor

For i In {0:#walls()-1}
  outSurf() -= {walls(i)}; 
  // Printf ("%g, ", walls(i));
EndFor

For i In {0:#outSurf()-1}
  Printf ("%g, ", outSurf(i));
EndFor

Printf("Number of wall surfaces = %g", #walls());
Printf("Number of farfield surfaces = %g", #outSurf());


// Walls
Physical Surface("Wall", 12) = walls();
Physical Surface("FarField", 20) = outSurf();
Physical Volume("Domain", 9) = vol();




// lf = 0.1; 
// lc = 4.0; 

l0 = 0.2; 
l1 = 0.3; 
l2 = 6.0; 

// lc = 6.0; 
// Meshing Related
// Give sizing to the wall points
// MeshSize{wpoints()} = lc;
// Mesh.MeshSizeFromCurvature = 20;

Field[1] = Distance;
Field[1].PointsList = {4, 18, 15};
Field[1].CurvesList = {10, 25};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = l0;
Field[2].SizeMax = l2;
Field[2].DistMin = 0.2;
Field[2].DistMax = 0.4;


Field[3] = Cylinder;
Field[3].Radius   = 2.5;
Field[3].VIn      = l1;
Field[3].VOut     = l2;
Field[3].XAxis    = 5;
Field[3].YAxis    = 0;
Field[3].ZAxis    = 0;
Field[3].XCenter  = 0.0;
Field[3].YCenter  =    0;
Field[3].ZCenter  =    0;


Field[4] = Ball;
Field[4].Radius   = 2.5;
Field[4].Thickness = 0.5;
Field[4].VIn      = l1;
Field[4].VOut     = l2;
Field[4].XCenter = -5.0;
Field[4].YCenter =    0;
Field[4].ZCenter =    0;


Field[5] = Ball;
Field[5].Radius   = 2.5;
Field[5].Thickness = 0.5;
Field[5].VIn      = l1;
Field[5].VOut     = l2;
Field[5].XCenter =  5.0;
Field[5].YCenter =    0;
Field[5].ZCenter =    0;


Field[6] = Ball;
Field[6].Radius   = 0.5;
Field[6].Thickness = 0.75;
Field[6].VIn      = l0;
Field[6].VOut     = l2;
Field[6].XCenter =  xmin;
Field[6].YCenter =    0;
Field[6].ZCenter =    0;


Field[7] = Min;
Field[7].FieldsList = {2, 3, 4, 5, 6};
Background Field = 7;





// // 2D mesh (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
// // 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of
// // Parallelograms, 11: Quasi-structured Quad)
// Mesh.Algorithm = 5; 

// // Mesh.Algorithm3D
// // 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9:
// // R-tree, 10: HXT)
// Mesh.Algorithm3D = 1; 





// // Activate optimizer
// Mesh.Optimize = 1; 
// // Optimize tetrahedra that have a quality below value
// Mesh.OptimizeThreshold = 0.3; 

// Mesh 2;

// Mesh 3;


// // Mesh output format (1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos
// // stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39:
// // inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab)
// Mesh.Format = 1; 
// Mesh.MshFileVersion = 2.2; 


// Save "test.msh";


 // // Mesh.....//
 // Mesh.MeshSizeFromCurvature = 10; // 10
 // Mesh.MeshSizeMin = 0.01;
 // Mesh.MeshSizeMax = 2.0;





























// For i In {0:#walls()-1}
//   MeshSize{ PointsOf {walls(i)} } = 0.15;  
// EndFor


// // Create Mesh
// If(MODEL==1 && BBOX==1)
//   MeshSize{ PointsOf{ Surface{2,3,4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18};}} = 0.15;  
//   MeshSize{ PointsOf{ Surface{1};}} = 5.0;  

//   Physical Surface("Wall", 12) = { 2,3,4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18};
//   Physical Surface("FarField", 20) = {1};
//   Physical Volume("Domain", 9) = {2};
// EndIf


// // Create Wall Mesh
// If(MODEL==1 && BBOX==2)
//   MeshSize{ PointsOf{ Surface{7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};}} = 0.15;  
//   MeshSize{ PointsOf{ Surface{1,2,3,4,5,6};}} = 5.0;  

//   Physical Surface("Wall", 12) = { 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
//   Physical Surface("FarField", 20) = {1,2,3,4,5,6};
//   Physical Volume("Domain", 9) = {2};
// EndIf


// // Create Wall Mesh
// If(MODEL==3 && BBOX==2)
//   MeshSize{ PointsOf{ Surface{7, 8, 9, 10, 11};}} = 0.05;  
//   MeshSize{ PointsOf{ Surface{1,2,3,4,5,6};}} = 5.0;  

//   Physical Surface("Wall", 12) = { 7,8,9,10,11};
//   Physical Surface("FarField", 20) = {1,2,3,4,5,6};
//   Physical Volume("Domain", 9) = {2};
// EndIf


// // // Mesh.....//
// Mesh.MeshSizeFromCurvature = 10; // 10
// Mesh.MeshSizeMin = 0.05;
// Mesh.MeshSizeMax = 3.0;

// Coherence;
// Mesh 2;
// Mesh 3;
// Mesh.Smoothing = 100;

// // // Mesh Visibility
// Mesh.SurfaceEdges = 1;
// Mesh.SurfaceFaces = 0;
// Mesh.VolumeEdges = 0;
