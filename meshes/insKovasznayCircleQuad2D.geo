res = DefineNumber[0.1];
Point(1) = {0, 0, 0, res};
Point(2) = {0, 1, 0, res};
Point(3) = {1, 0, 0, res};
Point(4) = {0, -1, 0, res};
Point(5) = {-1, 0, 0, res};

Circle(1) = {2, 1, 5};
Circle(2) = {5, 1, 4};
Circle(3) = {4, 1, 3};
Circle(4) = {3, 1, 2};

Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Line{1}; Line{2}; Line{3}; Line{4}; 
}
Line Loop(1) = {4, 1, 2, 3};

Plane Surface(1) = {1};
Physical Surface("Domain",9) = {1};
Physical Line("Inflow",2) = {1, 2, 4};
Physical Line("Outflow",3) = {3};
Recombine Surface {1};
