SetFactory("OpenCASCADE");
lc1  = DefineNumber[0.075];
lc2  = DefineNumber[5.00];

// scale_factor = 1.0/(63.5); 
scale_factor = 1/8.5; 

LoD = 1.0; 
xmin = 0.0; ymin = 0.0; zmin = 0.0;
xmax = 1.0; ymax = 1.0; zmax = 1.0;

// Printf("Loading model geometry with id = %g", MODEL);
v() = ShapeFromFile("./step_files/bodyalone.stp");
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
RI        = 50.0;
C0        = 10.0;  
Sphere(2) = {C0, 0, 0, RI, -Pi/2, Pi/2, 2*Pi};
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }


vol()     = Volume In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
walls()   = Surface In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};
outSurf() = Surface In BoundingBox {-1000, -1000, -1000, 1000, 1000, 1000};
wpoints() = Point In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};
wline()   = Curve In BoundingBox{xmin, ymin, zmin, xmax, ymax, zmax};

// exctract walls from all surfaces to get farfield
For i In {0:#walls()-1}
  outSurf() -= {walls(i)}; 
EndFor

// Volume ID check
For i In {0:#vol()-1}
  Printf ("volume Cehck : id of # volume = %g  over total of %g", vol(i),#vol() );
EndFor

// Point ID check
For i In {0:#wpoints()-1}
  Printf ("Point Cehck : id of # point = %g  over total of %g", wpoints(i),#wpoints() );
EndFor

// Curve ID check
For i In {0:#wline()-1}
  Printf ("Line Cehck : id of # line = %g  over total of %g", wline(i),#wline() );
EndFor


// Report Wall and Farfield geometric ids
For i In {0:#outSurf()-1}
  Printf ("FarField: %g", outSurf(i));
EndFor

For i In {0:#walls()-1}
  Printf ("Walls: %g", walls(i));
EndFor

Field[1] = Distance;
Field[1].PointsList = {1, 2, 3, 4, 5};
Field[1].CurvesList = {1,3,5,6,7,8,9,10};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc1;
Field[2].SizeMax = lc2;
Field[2].DistMin = 1.0;
Field[2].DistMax = RI;

Field[3] = Ball;
Field[3].Radius    =  0.15;
Field[3].Thickness =  1.0;
Field[3].VIn       =  0.1*lc1;
Field[3].VOut      =  lc2;
Field[3].XCenter   =  0.0;
Field[3].YCenter   =  0.0;
Field[3].ZCenter   =  0.0;

Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 4;

Mesh.Algorithm3D = 1; 
Mesh.Optimize = 1; 
