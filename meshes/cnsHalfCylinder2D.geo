Mesh.LcIntegrationPrecision = 1.e-2;
lc   = DefineNumber[0.02];
lc2  = DefineNumber[0.4];

R = DefineNumber[0.5]; 

xc   = DefineNumber[-0.5];
yc   = DefineNumber[0.0];

Point(1) = { xc+R, yc+R,  0.0, 0.5*lc}; 
Point(2) = { xc+R, yc,      0.0, lc}; 
Point(3) = { xc+R, yc-R,  0.0, 0.5*lc}; 

Circle(1) = {1, 2, 3};


Point(4) = { xc+R, 2.0,  0.0, lc}; 
Point(5) = { xc+R, -2.0,  0.0, lc}; 

Circle(2) = {4, 2, 5};
//+
Line(3) = {1, 4};
//+
Line(4) = {3, 5};
//+
Curve Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};
