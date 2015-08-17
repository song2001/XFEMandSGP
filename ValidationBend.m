% 2D bending of a beam
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=ABAQUSmesh2D;
%
% Material properties
% MAT:0(Linear elastic) / MAT:1 (Elastoplastic) / MAT:2 (CMSG plasticity)
% PROP=[YOUNG NU] (MAT=0)
% PROP=[YOUNG NU SYIELD] (MAT=1) 
% PROP=[YOUNG NU SYIELD ELLE N] (MAT=2)
MAT=0;
PROP=[260000, 0.3];
%
% Prescribed displacements [Node, DOF, Value]
Dnodes=[88 2 0;  % u2=0 in the middle node
175 1 0.9; 492 1 0.75; 150 1 0.6; 443 1 0.45; 125 1 0.3; 394 1 0.15; 100 1 0; 
% u1 in the upper right side
25 1 -0.9; 247 1 -0.75; 50 1 -0.6; 296 1 -0.45; 75 1 -0.3; 345 1 -0.15;
% u1 in the bottom right side
151 1 -0.9; 447 1 -0.75; 126 1 -0.6; 398 1 -0.45; 101 1 -0.3; 349 1 -0.15;
% u1 in the upper left side
1 1 0.9; 179 1 0.75; 26 1 0.6; 251 1 0.45; 51 1 0.3; 300 1 0.15; 76 1 0];
%
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[1, 1E-6, 30];    
% Calling main function 
EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');

