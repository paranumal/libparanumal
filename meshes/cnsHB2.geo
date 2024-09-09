lc   = DefineNumber[0.15];
lc2  = DefineNumber[0.5];
lc3  = DefineNumber[1.75];
fac  = DefineNumber[0.5];
// Note that Length = 8.5 x D for AGARD B
SetFactory("OpenCASCADE");

LoD = 1.0; 
xmin = 0.0; ymin = 0.0; zmin = 0.0;
xmax = 1.0; ymax = 1.0; zmax = 1.0;

// Printf("Loading model geometry with id = %g", MODEL);
v() = ShapeFromFile("./step_files/HB2.stp");
LoD = 8.5;
// Get the bounding box of the volume:
bbox() = BoundingBox Volume{v()};
xmin = bbox(0); ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3); ymax = bbox(4);zmax = bbox(5);
// Report the bounding box
For i In {0:#bbox()-1}
  // Printf("xmin = %g ", xmin);
  Printf ("before scaling bounding box x_%g= %g", i, bbox(i));
EndFor

//Set diameter to 1 in scaled geometry
// scale_factor = 1.0/(ymax - ymin); 
scale_factor = 1.0/(63.5); 

// Scale the volume so that L = 8.5; 
Dilate {{0, 0, 0}, {scale_factor, scale_factor, scale_factor}} { Volume{1};}

// Get the bounding box of the volume:
bbox() = BoundingBox Volume{v()};
xmin = bbox(0); ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3); ymax = bbox(4);zmax = bbox(5);
// Report the bounding box
For i In {0:#bbox()-1}
  // Printf("xmin = %g ", xmin);
  Printf ("before scaling bounding box x_%g= %g", i, bbox(i));
EndFor

// Get the bounding box of the volume:
bbox() = BoundingBox Volume{v()};
xmin = bbox(0); ymin = bbox(1);zmin = bbox(2);
xmax = bbox(3); ymax = bbox(4);zmax = bbox(5);

// Report the bounding box
For i In {0:#bbox()-1}
  Printf ("After Scaling: bounding box x_%g= %g", i, bbox(i));
EndFor

// RI = 3.0*(xmax - xmin)/2; 
RI        = 25.0;
C0        = -1.0;  
LC        = 75.0;  
Sphere(2) = {C0, 0, 0, RI, -Pi/2, Pi/2, 2*Pi};
Cylinder(3) = {C0, 0, 0, LC, 0, 0, RI, 2*Pi};
BooleanUnion{ Volume{3}; Delete; }{ Volume{2}; Delete; }
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }


vol()     = Volume In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
walls()   = Surface In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};
outSurf() = Surface In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
wpoints() = Point In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};

//Volume ID check
For i In {0:#vol()-1}
  Printf ("volume Cehck : id of # volume = %g  over total of %g", vol(i),#vol() );
EndFor

// exctract walls from all surfaces to get farfield
For i In {0:#walls()-1}
  outSurf() -= {walls(i)}; 
EndFor

For i In {0:#outSurf()-1}
  Printf ("FarField: %g", outSurf(i));
EndFor

For i In {0:#walls()-1}
  Printf ("Walls: %g", walls(i));
EndFor


// // Walls
// Physical Surface("Wall", 11) = walls();
// Physical Surface("FarField", 20) = outSurf();
// Physical Volume("Domain", 9) = vol();

Field[1] = Distance;
Field[1].PointsList = {2,1,4,5,6,3,2};
Field[1].CurvesList = {1,3,8,2,9,4,7,10,6,5,28,29};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc;
Field[2].SizeMax = lc3;
Field[2].DistMin = 1.0;
Field[2].DistMax = 2.0;

Field[3] = Cylinder;
Field[3].Radius   = 1.0;
Field[3].VIn      = fac*lc;
Field[3].VOut     = lc3;
Field[3].XAxis    = 3.0;
Field[3].YAxis    = 0.0;
Field[3].ZAxis    = 0.0;
Field[3].XCenter  = 3.0;
Field[3].YCenter  = 0.0;
Field[3].ZCenter  = 0.0;


Field[4] = Ball;
Field[4].Radius    =  1.25;
Field[4].Thickness =  0.75;
Field[4].VIn       =  fac*lc;
Field[4].VOut      =  lc3;
Field[4].XCenter   =  xmin+1.0;
Field[4].YCenter   =  0.0;
Field[4].ZCenter   =  0.0;

Field[5] = Ball;
Field[5].Radius    =  1.75;
Field[5].Thickness =  0.5;
Field[5].VIn       =  fac*lc;
Field[5].VOut      =  lc3;
Field[5].XCenter   =  xmax-0.5;
Field[5].YCenter   =  0.0;
Field[5].ZCenter   =  0.0;

Field[10] = Min;
Field[10].FieldsList = {2, 3, 4, 5,6};
Background Field = 10;

Mesh.Algorithm3D = 1; 
Mesh.Optimize = 1; 