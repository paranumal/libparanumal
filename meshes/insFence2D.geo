//+r1 = DefineNumber[0.0125];
//+r2 = DefineNumber[0.1];
//+r3 = DefineNumber[0.2];
//+r4 = DefineNumber[0.4];

r1 = DefineNumber[0.02];
r2 = DefineNumber[0.05];
r3 = DefineNumber[0.2];
r4 = DefineNumber[0.2];


//+
ffxmin = DefineNumber[-0.05];
ffxmax = DefineNumber[0.05];
ffymin = DefineNumber[0.0];
ffymax = DefineNumber[0.5];

//+
fbfacx = DefineNumber[8.0];
fbfacy = DefineNumber[2.5];

fbxmin = DefineNumber[fbfacx*ffxmin];
fbxmax = DefineNumber[fbfacx*ffxmax];
fbymin = DefineNumber[fbfacy*ffymin];
fbymax = DefineNumber[fbfacy*ffymax];

//+
fdxmin = DefineNumber[-2.5];
fdxmax = DefineNumber[5.0];
fdymin = DefineNumber[0.0];
fdymax = DefineNumber[2.5];

//+
Point(1) = {ffxmin, ffymin, 0.0, r1};
Point(2) = {ffxmin, ffymax, 0.0, r1};
Point(3) = {ffxmax, ffymax, 0.0, r1};
Point(4) = {ffxmax, ffymin, 0.0, r1};

//+
Point(5) = {fbxmin, fbymin, 0.0, r1};
Point(6) = {fbxmin, fbymax, 0.0, r1};
Point(7) = {fbxmax, fbymax, 0.0, r1};
Point(8) = {fbxmax, fbymin, 0.0, r1};

//+
Point(9)  = {fdxmin, fdymin, 0.0, r3};
Point(10) = {fdxmin, fdymax, 0.0, r4};
Point(11) = {fdxmax, fdymax, 0.0, r4};
Point(12) = {fdxmax, fdymin, 0.0, r2};

//+
Point(13) = {fdxmin, fbymax, 0.0, r3};
Point(14) = {fdxmax, fbymax, 0.0, r2};
Point(15) = {fbxmin, fdymax, 0.0, r4};
Point(16) = {fbxmax, fdymax, 0.0, r4};

//+
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 4};
Line(9) = {8, 12};
Line(10) = {12, 14};
Line(11) = {14, 7};
Line(12) = {11, 14};
Line(13) = {16, 7};
Line(14) = {16, 11};
Line(15) = {16, 15};
Line(16) = {15, 6};
Line(17) = {6, 13};
Line(18) = {13, 10};
Line(19) = {10, 15};
Line(20) = {13, 9};
Line(21) = {9, 5};
Line Loop(22) = {5, 6, 7, 8, 1, 2, 3, 4};
Plane Surface(23) = {22};
Line Loop(24) = {7, 9, 10, 11};
Plane Surface(25) = {24};
Line Loop(26) = {11, -13, 14, 12};
Plane Surface(27) = {26};
Line Loop(28) = {16, 6, -13, 15};
Plane Surface(29) = {28};
Line Loop(30) = {16, 17, 18, 19};
Plane Surface(31) = {30};
Line Loop(32) = {17, 20, 21, 5};
Plane Surface(33) = {32};

Physical Line("Wall",1) = {3, 2, 1, 21, 4, 8, 9};
Physical Line("Inflow",2) = {20, 18};
Physical Line("Outflow",3) = {12, 10};
Physical Line("YSlip",5) = {19, 15, 14};
Physical Surface("Domain",9) = {33, 31, 23, 29, 27, 25};
