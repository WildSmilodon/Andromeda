//
//
// number of layers
nx=50;
ny=10;
nz=20;
s = 0.2;
// length
lx=30;
ly=6;
lz=2;

Point(1) = {0,-ly/2,-lz/2,s};
Point(2) = {0,+ly/2,-lz/2,s};
Point(3) = {0,+ly/2,0,s};
Point(4) = {0,-ly/2,0,s};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};

Point(5) = {0,-ly/2,+lz/2,s};
Point(6) = {0,+ly/2,+lz/2,s};

Line(5) = {3, 6};
Line(6) = {6, 5};
Line(7) = {5, 4};

Line Loop(2) = {-3, 5, 6, 7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// outlet side

Point(11) = {lx,-ly/2,-lz/2,s};
Point(12) = {lx,+ly/2,-lz/2,s};
Point(13) = {lx,+ly/2,0,s};
Point(14) = {lx,-ly/2,0,s};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};
Line Loop(11) = {11, 12, 13, 14};

Point(15) = {lx,-ly/2,+lz/2,s};
Point(16) = {lx,+ly/2,+lz/2,s};

Line(15) = {13, 16};
Line(16) = {16, 15};
Line(17) = {15, 14};

Line Loop(12) = {-13, 15, 16, 17};

Plane Surface(11) = {11};
Plane Surface(12) = {12};

Line(21) = {6,16};
Line(22) = {3,13};
Line(23) = {2,12};
Line(24) = {1,11};
Line(25) = {4,14};
Line(26) = {5,15};

Line Loop (21) = { -23, 2, 22, -12 };
Line Loop (22) = { -22, 5, 21, -15 };
Line Loop (23) = { 24, -14, -25, 4 };
Line Loop (24) = { 25, -17, -26, 7 };

Plane Surface(21) = {21};
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};


Line Loop (31) = { 26, -16, -21, 6 };
Line Loop (32) = { 24, 11, -23, -1 };

Plane Surface(31) = {31};
Plane Surface(32) = {32};


Transfinite Line {21,22,23,24,25,26} = nx+1 Using Progression 1.05;
Transfinite Line {1, 3, 6} = ny+1 Using Progression 1;
Transfinite Line {2, 4, 5, 7} = nz+1 Using Bump 0.25;
Transfinite Line {11, 13, 16} = ny+1 Using Progression 1;
Transfinite Line {12, 14, 15, 17} = nz+1 Using Bump 0.25;

Recombine Surface {1};
Transfinite Surface {1};
Recombine Surface {2};
Transfinite Surface {2};

Recombine Surface {11};
Transfinite Surface {11};
Recombine Surface {12};
Transfinite Surface {12};

Recombine Surface {21};
Transfinite Surface {21};
Recombine Surface {22};
Transfinite Surface {22};

Recombine Surface {23};
Transfinite Surface {23};
Recombine Surface {24};
Transfinite Surface {24};

Recombine Surface {31};
Transfinite Surface {31};
Recombine Surface {32};
Transfinite Surface {32};

// CASE : NORMALS OUT OF DOMAIN
// Reverse normals
ReverseMesh Surface{1}; // this was done using "-" in Physical Sruface in v2 mesh, in v4 it is not
ReverseMesh Surface{2}; // this was automatic in v2 mesh, in v4 it is not
ReverseMesh Surface{32}; // this was automatic in v2 mesh, in v4 it is not


Physical Surface("inlet") = {2};
Physical Surface("step") = {1};
Physical Surface("outlet") = {11,12};
Physical Surface("top") = {31};
Physical Surface("bottom") = {32};
Physical Surface("front") = {23,24};
Physical Surface("back") = {21,22};//+
