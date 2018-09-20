// coarse = DefineNumber[12.5];
// fine   = DefineNumber[0.625];
// fineC  = DefineNumber[0.10];

coarse = DefineNumber[1.0];
fine   = DefineNumber[0.2];
fineC  = DefineNumber[0.10];

xmax  = DefineNumber[25];
xmin  = DefineNumber[-16];
ymax  = DefineNumber[22];
ymin  = DefineNumber[-22];

xcmax  = DefineNumber[0.5];
xcmin  = DefineNumber[-0.5];
ycmax  = DefineNumber[0.5];
ycmin  = DefineNumber[-0.5];

xbmax  = xmax;
xbmin  = DefineNumber[-2];
ybmax  = DefineNumber[2];
ybmin  = DefineNumber[-2];


Point(1) = {xmin, ymin, 0, coarse};
Point(2) = {xmin, ymax, 0, coarse};
Point(3) = {xmax, ymax, 0, coarse};
Point(4) = {xmax, ymin, 0, coarse};

Point(6) = {xcmin, ycmax, 0, fineC};
Point(7) = {xcmax, ycmax, 0, fineC};
Point(8) = {xcmin, ycmin, 0, fineC};
Point(9) = {xcmax, ycmin, 0, fineC};

Point(10) = {xbmin,ybmin, 0, fine};
Point(11) = {xbmin, ybmax, 0, fine};
Point(12) = {xbmax, ybmax, 0, fine};
Point(13) = {xbmax,ybmin, 0, fine};

Line(1) = {2, 3};
Line(2) = {12, 3};
Line(3) = {11, 12};
Line(4) = {2, 1};
Line(5) = {4, 1};
Line(6) = {13, 4};
Line(7) = {13, 10};
Line(8) = {11, 10};
Line(9) = {12, 13};
Line(10) = {6, 7};
Line(11) = {7, 9};
Line(12) = {9, 8};
Line(13) = {8, 6};

Line Loop(14) = {1, -2, -3, 8, -7, 6, 5, -4};
Plane Surface(15) = {14};
Line Loop(16) = {3, 9, 7, -8};
Line Loop(17) = {10, 11, 12, 13};
Plane Surface(18) = {16, 17};

Coherence;
Physical Line("Wall",1) = {12, 11, 13, 10};
Physical Line("Inflow",2) = {4,1,5};
Physical Line("Outflow",3) = {2, 9, 6};
Physical Surface("Domain",9) = {15, 18};
Coherence;
Coherence;
