// Gmsh project created on Mon Jul 27 10:39:02 2020
SetFactory("OpenCASCADE");

// Import processed cell surface
//Merge "larger_sphere.msh";
//ClassifySurfaces{Pi/4, 0, 0};
//CreateGeometry;
//Surface Loop(43)={0};

// Resolution control
tol = .05;
avg_len = 4;

lc = 10;

// Define dimensions of imaging box
xy_bound = 149.5;
z_bound = 120;

cx = 74.5;
cy = 74.5;
cz = 60;
rd = 30;


//Point 2
P2x = cx + rd;
P2y = cy;
P2z = cz;

//Point 3
P3x = cx;
P3y = cy + rd;
P3z = cz;

//Point 4
P4x = cx - rd;
P4y = cy;
P4z = cz;

//Point 5
P5x = cx;
P5y = cy - rd;
P5z = cz;


//Point 6
P6x = cx;
P6y = cy;
P6z = cz + rd;

//Point 7
P7x = cx;
P7y = cy;
P7z = cz - rd;

// Create Sphere manually


// Construct imaging box with GMSH's basic geometry objects
Point(1) = {0, 0, 0, lc};
Point(2) = {0, xy_bound, 0, lc};
Point(3) = {xy_bound, xy_bound, 0, lc};
Point(4) = {xy_bound, 0, 0, lc};
Point(5) = {0, 0, z_bound, lc};
Point(6) = {0, xy_bound, z_bound, lc};
Point(7) = {xy_bound, xy_bound, z_bound, lc}; 
Point(8) = {xy_bound, 0, z_bound, lc};
// Points for Circle
Point(9) = {cx, cy, cz, lc};
Point(10) = {P2x, P2y, P2z, lc};
Point(11) = {P3x, P3y, P3z, lc};
Point(12) = {P4x, P4y, P4z, lc};
Point(13) = {P5x, P5y, P5z, lc};
Point(14) = {P6x, P6y, P6z, lc};
Point(15) = {P7x, P7y, P7z, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,5};
Line(10) = {2, 6};
Line(11) = {3,7};
Line(12) = {4,8};
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {1, 10, -6, -5};
Curve Loop(3) = {2, 11, -7, -10};
Curve Loop(4) = {3, 12, -8, -11};
Curve Loop(5) = {6, 7, 8, 9};
Curve Loop(6) = {5, -9, -12, 4};

// Circle in XY Plane
Circle(13) = {10, 9, 11};
Circle(14) = {11, 9, 12};
Circle(15) = {12, 9, 13};
Circle(16) = {13, 9, 10};

// Circle in XZ Plane
Circle(17) = {14, 9, 10};
Circle(18) = {10, 9, 15};
Circle(19) = {15, 9, 12};
Circle(20) = {12, 9, 14};

// Circle in YZ Plane
Circle(21) = {14, 9, 13};
Circle(22) = {13, 9, 15};
Circle(23) = {15, 9, 11};
Circle(24) = {11, 9, 14};

// Create Sphere Surfaces / Create Shell
Line Loop(13) = {13, 24, 17};
Surface(14) = {13};

Line Loop(15) = {24, -20, -14};
Surface(16) = {15};

Line Loop(17) = {20, 21, -15};
Surface(18) = {17};

Line Loop(19) = {21, 16, -17};
Surface(20) = {19};

Line Loop(21) = {13, 23, -18};
Surface(22) = {21};

Line Loop(23) = {22, -18, -16};
Surface(24) = {23};

Line Loop(25) = {22, 19, 15};
Surface(26) = {25};

Line Loop(27) = {23, 14, -19};
Surface(28) = {27};

// Create a Solid
Surface Loop(29) = {16, 14, 22, 28, 26, 24, 20, 18};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Surface Loop(42) = {1:6}; // Ultimate imaging box object


Volume(100) = {42, 29}; // Gel mesh


// Create sphere
//Sphere(1) = {74.5, 74.5, 60, 30, -Pi/2, Pi/2, 2*Pi};
//Plane Surface(8) = {13};
//Surface Loop(44) = {8};



// Resolution
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0.05;

// definitely has an effect
Mesh.CharacteristicLengthMax = avg_len + tol;
Mesh.CharacteristicLengthMin = avg_len - tol;
Mesh.CharacteristicLengthFactor = 1;

