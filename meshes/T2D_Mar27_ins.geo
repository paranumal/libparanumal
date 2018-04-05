
lc = 1.00;

// x-bounds
xmin  = -5.0;
xmax  =  13.0;

// y-bounds
ymin  =  -6.0;
ymax  =  6.0;

// z-bounds
zall = 0.0;

h = 1.0;            // height of fence
w = 0.1;            // width of fence
t = 5.0;            //legth of tail

fx1 = 0.00;         // front of fence
fx2 = fx1 + w;      // width of fence
fx3 = fx2 + t;      // legth of tail

fy1 = -w/2;         // base of tail
fy2 =  w/2;         // top of tail
fy3 =  fy1-h;       // base of fence
fy4 =  fy2+h;       // base of fence

//------ Domain-------------------------
Point(1) = { xmin,  ymin, zall, lc*1.0};
Point(2) = { xmax,  ymin, zall, lc*1.0};
Point(3) = { xmin,  ymax, zall, lc*1.0};
Point(4) = { xmax,  ymax, zall, lc*1.0};

//---------------------------------------
Point(5) = { fx1,  fy3, zall, 0.025};  // base front
Point(6) = { fx1,  fy4, zall, 0.025};  // top  front

Point(7) = { fx2,  fy4, zall, 0.025};  // top rear 
Point(8) = { fx2,  fy2, zall, 0.050};  // top rear base

Point(9)  = { fx3,  fy2, zall, 0.050};  // tail top 
Point(10) = { fx3,  fy1, zall, 0.050};  // tail base

Point(11) = { fx2,  fy1, zall, 0.025};  // bottom rear base
Point(12) = { fx2,  fy3, zall, 0.050};  // bottom rear


//---------------------------------------
//Domain
Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {4,3};
Line(4) = {3,1};

//Fence
Line(5)  = { 5, 6};
Line(6)  = { 6, 7};
Line(7)  = { 7, 8};
Line(8)  = { 8, 9};
Line(9)  = { 9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12, 5};


//---------------------------------------
// main fluid domain
//---------------------------------------
LLDomain=newll; Line Loop(LLDomain) = {1,2,3,4};
LLFence=newll;  Line Loop(LLFence) = {5,6,7,8,9,10,11,12};
PSDomain=news;  Plane Surface(PSDomain) = {LLDomain, LLFence};

Physical Surface("Domian",9) = {PSDomain};
Physical    Line("Wall", 1) = {5,6,7,8,9,10,11,12};
Physical    Line("Inflow",2) = {4};
Physical    Line("Outflow",3) = {2};
Physical    Line("YSlip",5) = {1,3};

//---------------------------------------------------------
// use Fields to tune wake density
//---------------------------------------------------------

Field[1] = Box;
Field[1].VIn  =  0.035;
Field[1].VOut =  0.50;
Field[1].XMin = -0.8;
Field[1].XMax =  2.5;
Field[1].YMin =  -1.8;
Field[1].YMax =   1.8;

Field[2] = Box;
Field[2].VIn  =  0.05;
Field[2].VOut =  0.50;
Field[2].XMin = -2.2;
Field[2].XMax =  5.5;
Field[2].YMin =  -2.7;
Field[2].YMax =   2.7;

Field[3] = Box;
Field[3].VIn  =  0.05;
Field[3].VOut =  0.50;
Field[3].XMin = -3.3;
Field[3].XMax =  0.0;
Field[3].YMin =  -1.0;
Field[3].YMax =   1.0;

//Field[4] = Box;
//Field[4].VIn  =  0.10;
//Field[4].VOut =  0.50;
//Field[4].XMin = -4.5;
//Field[4].XMax =  8.0;
//Field[4].YMin =  0.0;
//Field[4].YMax =  3.0;

Field[5] = Box;
Field[5].VIn  =  0.15;
Field[5].VOut =  0.50;
Field[5].XMin = -4.0;
Field[5].XMax = 12.0;
Field[5].YMin = -4.0;
Field[5].YMax =  4.0;

Field[6] = Box;
Field[6].VIn  =  0.02;
Field[6].VOut =  0.50;
Field[6].XMin = -0.2;
Field[6].XMax =  0.5;
Field[6].YMin =  0.8;
Field[6].YMax =  1.2;

Field[7] = Box;
Field[7].VIn  =  0.02;
Field[7].VOut =  0.50;
Field[7].XMin = -0.2;
Field[7].XMax =  0.5;
Field[7].YMin =  -1.2;
Field[7].YMax =  -0.8;


Field[8] = Min;
Field[8].FieldsList = {1,2,3,4,5,6,7};
Background Field = 8;


