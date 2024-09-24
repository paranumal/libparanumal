Mesh.LcIntegrationPrecision = 1.e-2;
lc   = DefineNumber[0.02];
lc2  = DefineNumber[0.4];

R = DefineNumber[0.5]; 

xc   = DefineNumber[0];
yc   = DefineNumber[0];

Point(1) = { xc, yc,  0.0, lc}; 
Point(2) = { xc+R, yc,  0.0, lc}; 
Point(3) = { xc-R, yc,  0.0, lc}; 

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};

xa   = DefineNumber[0.00];
ya   = DefineNumber[10.0];

xb   = DefineNumber[-4.0];
yb   = DefineNumber[0.0];

xc   = DefineNumber[0.00];
yc   = DefineNumber[-10.0];

xd   = DefineNumber[10.0];
yd   = DefineNumber[0.0];

Point(80) = { xc, yc,  0.0, lc}; 
Point(81) = { xa, ya,  0.0, lc2}; 
Point(82) = { xb, yb,  0.0, lc2}; 
Point(83) = { xc, yc,  0.0, lc2}; 
Point(84) = { xd, yd,  0.0, lc2}; 

Circle(3) = {81, 82, 83};
Circle(4) = {83, 84, 81};

// //+
Curve Loop(1)    = {1,2};
Curve Loop(2)    = {3, 4};
Plane Surface(1) = {1, 2};

// Physical Line("wall",1)      = {1,2};
// Physical Line("far",2)       = {3,4};
// Physical Surface("fluid",3)  = {1};