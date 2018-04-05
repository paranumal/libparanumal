// K=18,832 [quads]

lc = 1.00;

// x-bounds
xmin  = -8.0;
xmax  =  13.0;


// y-bounds
ymin  =  0.0;
ymax  =  6.0;

// z-bounds
zall = 0.0;

fx1 = 0.00;         // front of fence
fx2 = fx1 + 0.10;   // width of fence
fy1 = 0.00;         // base of fence
fy2 = 1.00;         // top of fence


Point(1) = { xmin,   ymin, zall, lc*1.0};
//---------------------------------------
Point(3) = { fx1,  fy1, zall, 0.050};  // base front
Point(4) = { fx1,  fy2, zall, 0.025};  // top  front
Point(5) = { fx2,  fy2, zall, 0.025};  // top  rear
Point(6) = { fx2,  fy1, zall, 0.050};  // base rear
//---------------------------------------
Point(8) = {xmax,   ymin,  zall, lc*1.0};

Point(10)= {xmax,   ymax,  zall, lc*1.0};
Point(13)= {xmin,   ymax,  zall, lc*1.0};
//---------------------------------------
Line(2) = {1,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,8};

Line(17)= {13, 1};
Line(19)= { 8,10};

Line(22)= {10,13};



//---------------------------------------
// main fluid domain
//---------------------------------------
LLDomain=newll; Line Loop(LLDomain) = {2,3,4,5,6,19,22,17};
PSDomain=news;  Plane Surface(PSDomain) = {LLDomain};

Physical Surface("Domian",9) = {PSDomain};
Physical    Line("Wall", 1) = {2,3,4,5,6};
Physical    Line("Inflow",2) = {17};
Physical    Line("Outflow",3) = {19};
Physical    Line("YSlip",5) = {22};

//---------------------------------------------------------
// use Fields to tune wake density
//---------------------------------------------------------

Field[1] = Box;
Field[1].VIn  =  0.035;
Field[1].VOut =  0.50;
Field[1].XMin = -0.8;
Field[1].XMax =  2.5;
Field[1].YMin =  0.0;
Field[1].YMax =  1.8;

Field[2] = Box;
Field[2].VIn  =  0.05;
Field[2].VOut =  0.50;
Field[2].XMin = -2.2;
Field[2].XMax =  5.0;
Field[2].YMin =  0.0;
Field[2].YMax =  2.7;

Field[3] = Box;
Field[3].VIn  =  0.05;
Field[3].VOut =  0.50;
Field[3].XMin = -3.3;
Field[3].XMax =  0.0;
Field[3].YMin =  0.0;
Field[3].YMax =  1.0;

Field[4] = Box;
Field[4].VIn  =  0.10;
Field[4].VOut =  0.50;
Field[4].XMin = -4.5;
Field[4].XMax =  8.0;
Field[4].YMin =  0.0;
Field[4].YMax =  3.0;

Field[5] = Box;
Field[5].VIn  =  0.15;
Field[5].VOut =  0.50;
Field[5].XMin = -6.0;
Field[5].XMax = 12.0;
Field[5].YMin =  0.0;
Field[5].YMax =  4.0;

Field[6] = Box;
Field[6].VIn  =  0.02;
Field[6].VOut =  0.50;
Field[6].XMin = -0.2;
Field[6].XMax =  0.5;
Field[6].YMin =  0.8;
Field[6].YMax =  1.2;


Field[7] = Min;
Field[7].FieldsList = {1,2,3,4,5,6};
Background Field = 7;
Recombine Surface {24};
