                                                                  
R=0.5;
R1=R/2;
R2=R/4;

kot=10;
kotRad=kot*Pi/180;

xC=R*Cos(kotRad);
yC=R*Sin(kotRad);
p=0.25;
Point(0) = { 0, 0, 0, p }; 

Point(1) = { -R, 0, 0,p }; 
Point(2) = { -xC, yC, 0,p }; 
Point(3) = { xC, yC, 0,p }; 
Point(4) = { R, 0, 0 ,p}; 

Point(5) = { -R1, 0, 0, p};
Point(6) = { 0, R1, 0, p }; 
Point(7) = { R1, 0, 0, p };

Point(8) = { -R2, 0, 0, p};
Point(9) = { 0, R2, 0, p }; 
Point(10) = { R2, 0, 0, p };


//+
Circle(0) = {6, 0, 5};
//+
Circle(1) = {7, 0, 6};
//+
Circle(2) = {1, 0, 2};
//+
Circle(3) = {2, 0, 3};
//+
Circle(4) = {3, 0, 4};

//+
Circle(5) = {9, 0, 8};
//+
Circle(6) = {10, 0, 9};

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{2}; Line{4}; Line{3};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{0}; Line{1};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{5}; Line{6};
}



Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{2}; Line{4}; Line{3};
}

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{0}; Line{1};
}

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{5}; Line{6};
}


Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{7}; Line{13}; Line{10};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{17}; Line{20};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{23}; Line{26};
}


Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
   Line{51}; Line{58}; Line{54};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{64}; Line{61};
}

Extrude {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Line{70}; Line{67};
}


Physical Surface("plus") = {12,-34,60,78};
//+
Physical Surface("minus") = {9,-31,53,75};
//+
Physical Surface("outer-wall") = {16,-38,57,82};

Physical Surface("middle-wall") = {-19,-22,44,41,-63,-66,-88,-85};


Physical Surface("inner-wall") = {50,-69,-72,-25,-28,47,-94,-91};

Mesh.CharacteristicLengthFactor = 0.5;
