//
//  gmsh mreza za kroglo v kanalu
//

// Novi gmsh ne zasuka normal, uporabi ali
// gmsh -2 -format msh2 sphere-in-channel.geo
// ali
// ukaz ReverseMesh Surface{16}; to naredi

// presek kanala
D = 2;

// Dolzina kanala
L = 6;

// Velikost elementa v kanalu
lcar=0.5;


s=0.1; // Velikost elementa v krogli
Mesh.CharacteristicLengthFromCurvature = 0.05;
r=0.5; // radij krogle



// *********************************************************************
// **************** KROGLA *********************************************
// *********************************************************************

Point(1) = {0, 0, 0, s};
Point(2) = {-r, 0, 0, s};
Point(3) = {r, 0, 0, s};
Point(4) = {0, -r, 0, s};
Point(5) = {0, r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};
Circle(1) = {3, 1, 5};
Circle(2) = {2, 1, 5};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 3};
Circle(5) = {6, 1, 3};
Circle(6) = {6, 1, 2};
Circle(7) = {2, 1, 7};
Circle(8) = {7, 1, 3};
Circle(9) = {5, 1, 6};
Circle(10) = {6, 1, 4};
Circle(11) = {4, 1, 7};
Circle(12) = {7, 1, 5};
Line Loop(13) = {9, 5, 1};
Surface(14) = {13};
Line Loop(15) = {2, 9, 6};
Surface(16) = {15};
Line Loop(17) = {7, 12, -2};
Surface(18) = {17};
Line Loop(19) = {12, -1, -8};
Surface(20) = {19};
Line Loop(21) = {6, 3, -10};
Surface(22) = {21};
Line Loop(23) = {10, 4, -5};
Surface(24) = {23};
Line Loop(25) = {11, -7, 3};
Surface(26) = {25};
Line Loop(27) = {8, -4, 11};
Surface(28) = {27};

// Reverse normals
ReverseMesh Surface{16}; // this was done using "-" in Physical Sruface in v2 mesh, in v4 it is not
ReverseMesh Surface{18}; // this was automatic in v2 mesh, in v4 it is not
ReverseMesh Surface{26}; // this was automatic in v2 mesh, in v4 it is not

//Physical surfaces of domain
Physical Surface("krogla") = {14, 16, 18, 20, 22, 24, 26, 28};


// *********************************************************************
// **************** KANAL ***********************************************
// *********************************************************************

//+
SetFactory("OpenCASCADE");

Box(1) = {-L/2, -D/2, -D/2, L, D, D};

//+
Physical Surface("inlet") = {29};
//+
Physical Surface("outlet") = {30};
//+
Physical Surface("front") = {31};
//+
Physical Surface("back") = {32};
//+
Physical Surface("bottom") = {33};
//+
Physical Surface("top") = {34};
