//Inputs for Gmsh
boxdim = 1;
gridsize = boxdim/2;

Point(1) = {0,0,0,gridsize};
Point(2) = {boxdim,0,0,gridsize};
Point(3) = {boxdim,boxdim,0,gridsize};
Point(4) = {0,boxdim,0,gridsize};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};

Plane Surface(10) = 9;

// Square structured mesh
Transfinite Line{5,6,7,8} = boxdim/gridsize;
Transfinite Surface{10};
Recombine Surface{10};
