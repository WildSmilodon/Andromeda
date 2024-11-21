//
//  gmsh mreza za valj ki se zozi v valj
//
SetFactory("OpenCASCADE");

// radij valja
r=0.5;
r2=0.5/3;
// dolzina conea
lcone=0.085714285714286;

// dolzina valja
l=2.142857142857143;
l2=0.857142857142857;
z2=l+lcone;



// Velikost elementa
lcar=0.035;

Point(1) = {0,0,0,lcar};
Point(2) = {r,0,0,lcar};
Point(3) = {0,r,0,lcar};
Point(4) = {-r,0,0,lcar};
Point(5) = {0,-r,0,lcar};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

//
Point(111) = {0,0,z2,lcar};
Point(112) = {r2,0,z2,lcar};
Point(113) = {0,r2,z2,lcar};
Point(114) = {-r2,0,z2,lcar};
Point(115) = {0,-r2,z2,lcar};
//
Circle(111) = {112,111,113};
Circle(112) = {113,111,114};
Circle(113) = {114,111,115};
Circle(114) = {115,111,112};

Line Loop(115) = {111,112,113,114};
Plane Surface(116) = {115};

Extrude {0,0,l2} {
    Surface{116};
}


Extrude {0,0,l} {
    Surface{6};
}
  

Cone(1111) = {0, 0, l, 0, 0, lcone, r, r2, 2*Pi};

ReverseMesh Surface{6}; // this was automatic in v2 mesh, in v4 it is not

Physical Surface("spredaj") = {6}; 
Physical Surface("stena") = {122,123,124,125,127,117,118,119,120};
Physical Surface("zadaj") = {121};
