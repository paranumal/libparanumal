Mesh.LcIntegrationPrecision = 1.e-2;
lc   = DefineNumber[0.025];
lc2  = DefineNumber[0.3];
lc3  = DefineNumber[1.0];

R = DefineNumber[0.5]; 


SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, R, 0, 2*Pi};


xc   = DefineNumber[0];
yc   = DefineNumber[0];

xi   = DefineNumber[0];
yi   = DefineNumber[2.0];

yomin   = DefineNumber[-5.0];
yomax   = -yomin; 
xomin   = DefineNumber[5.0];
// xomax   = DefineNumber[5.0];
xomax   = DefineNumber[5.0];
// External box
Point(800) = { xc, yc,  0.0, lc3}; 

Point(801) = { xc, yomin,  0.0, lc3}; 
Point(802) = { xomax, yomin,  0.0, lc3}; 
Point(803) = { xomax, yomax,  0.0, lc3}; 
Point(804) = { xc,   yomax,  0.0, lc3}; 
Point(805) = { yomin,  yc,  0.0, lc3}; 

// Internal region
Point(806) = { xi, -yi,  0.0, lc2}; 
Point(807) = { xi,  yi,  0.0, lc2}; 
Point(808) = { -yi,  0.0,  0.0, lc2}; 

Point(809) = { xomax, -(yi+0.5),  0.0, lc2}; 
Point(810) = { xomax,  (yi+0.5),  0.0, lc2}; 


Line(5)  = {801,802}; 
Line(6)  = {802,809}; 
Line(7)  = {809,810}; 
Line(8)  = {810,803}; 
Line(9)  = {803,804}; 

Circle(10) = {804, 800, 805};
Circle(11) = {805, 800, 801};


Line(12)  = {806,809}; 
Line(13)  = {807,810}; 

Circle(14) = {807, 800, 808};
Circle(15) = {808, 800, 806};


Curve Loop(1) = {1};
Curve Loop(2) = {15, 12, 7, -13, 14};
Curve Loop(3) = {11, 5, 6, 7, 8, 9, 10};


Plane Surface(1) = {2, 3};
Plane Surface(2) = {1, 2};



Field[1] = Distance;
Field[1].EdgesList = {1};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc;
Field[2].SizeMax = lc2;
Field[2].DistMin = 0.2;
Field[2].DistMax = 0.3;

// Let's use the minimum of all the fields as the background mesh size field
Field[4] = Min;
Field[4].FieldsList = {1,2};
Background Field = 2;

Coherence;

// We are using Geometric Identites NOW!!!!

// Isothermall Wall
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
//+

