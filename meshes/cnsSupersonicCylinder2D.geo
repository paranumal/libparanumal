Mesh.LcIntegrationPrecision = 1.e-2;
lc   = DefineNumber[0.020];
lc2  = DefineNumber[0.40];
lc3  = DefineNumber[1.5];

R = DefineNumber[0.5]; 

SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, R, 0, 2*Pi};

xc   = DefineNumber[0];
yc   = DefineNumber[0];

xa   = DefineNumber[0.00];
ya   = DefineNumber[10.0];

xb   = DefineNumber[-4.0];
yb   = DefineNumber[0.0];

xc   = DefineNumber[0.00];
yc   = DefineNumber[-10.0];

xd   = DefineNumber[10.0];
yd   = DefineNumber[0.0];

Point(80) = { xc, yc,  0.0, lc}; 
Point(81) = { xa, ya,  0.0, lc3}; 
Point(82) = { xb, yb,  0.0, lc3}; 
Point(83) = { xc, yc,  0.0, lc3}; 
Point(84) = { xd, yd,  0.0, lc3}; 

Circle(2) = {81, 82, 83};
Circle(3) = {83, 84, 81};

//+
Curve Loop(1)    = {1};
Curve Loop(2)    = {2, 3};
Plane Surface(1) = {1, 2};



Field[1] = Distance;
Field[1].EdgesList = {1};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc;
Field[2].SizeMax = lc2;
Field[2].DistMin = 0.4;
Field[2].DistMax = 0.5;

// Let's use the minimum of all the fields as the background mesh size field
Field[4] = Min;
Field[4].FieldsList = {1,2};
Background Field = 2;
Coherence;

// We are using Geometric Identites NOW!!!!

// // Isothermall Wall
// Physical Line("Wall",12)      = {1, 2, 3, 4};
// Physical Line("Drichlet",41)  = {5,6, 7, 8, 9,10,11};
// Physical Surface("Domain",9)  = {1,2};//+

// // Use geometric IDs and let user change it to BC IDs
// Physical Line("Wall_1",1)      = {1};
// Physical Line("Wall_2",2)      = {2};
// Physical Line("Wall_3",3)      = {3};
// Physical Line("Wall_4",4)      = {4};
// Physical Line("Drichlet_1",5) = {10,11};
// Physical Line("Drichlet_2",6)   = {5,9};
// Physical Line("Drichlet_3",7)   = {6, 7, 8};
// Physical Surface("Domain",9)  = {1,2};//+
// +

