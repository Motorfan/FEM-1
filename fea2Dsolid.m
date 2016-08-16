function fea2Dsolid
%FEA 2D code for plane stress by Alex P. Martinez
%Meshing function modified from MESHDEMO2d Distmesh2d examples.
%Other functions used with MESHDEMO2d are unmodified.
%Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

%Generate the same random numbers to get same results
rand('state',111);
set(gcf,'rend','z');

E         = 200e9; %In Pascals
nu        = 0.3;   %Dimensionless
%Applied distributed load per unit thickness 
sigma     = 1e4;   %Has same units as pressure (Pascals) 
r         = 0.05;  %Radius of hole (~0.05-0.6; or distmesh doesn't work)


%XXXXXXXXXXXXXXXXXXXXXXX     PREPROCESSOR         XXXXXXXXXXXXXXXXXXXXXXXXX

%GENERATE MESH
disp(' ')
disp('FINITE ELEMENT ANALYSIS OF A HOLE IN A PLATE')

disp(' ')
disp('Generating mesh')

%Distance function
str = sprintf('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,%d))', r);
fd=inline(str,'p');
%Relative desired edge length function
fh=inline('6*sqrt(sum(p.^2,2))-0.2','p');
%Distance between points in initial distribution
h0=0.004;
%Bounding box: to bound whole shape!
box=[-1,-1;1,1];
%Positions of fixed nodes (for meshing)
fix=[-1,-1;-1,1;1,-1;1,1];
%Main mesh generating function
[p,t]=distmesh2d(fd,fh,h0,box,fix);
%Remove duplicated/unused nodes and place element nodes in 
%counterclockwise order to get nonsingular jacobians
[p,t] = fixmesh(p,t);
%Checking quality
q=simpqual(p,t);
u=uniformity(p,t,fh);
%Print mesh analysis
disp(sprintf(' - Mesh min quality %.2f',   min(q)))
disp(sprintf(' - Mesh uniformity  %.1f%%',  100*u))

elemType  = 'T3';
gaussType = 'TRIANGULAR';
gaussDeg  = 1;
print     = 1;
%Preprocessor and processor
[K, U] = fep(E, nu, sigma, r, p, t, elemType, gaussType, gaussDeg, print);
%Postprocessor
[Sx, Sy, Sxy] = fepp(E, nu, p, t, U, elemType, gaussType, print);
%Graph displacements
graphDs(elemType, p, t, K, U);
%Graph stresses
%graphSs(elemType, p, t, Sx );
graphSs(elemType, p, t, Sy );
%graphSs(elemType, p, t, Sxy);