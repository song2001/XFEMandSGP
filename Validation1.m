% 2D multiple elements uniaxial displacement controlled test
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
% Introduce manually the number of nodes per element:
Q=4;
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=ABAQUSmesh2D(Q);
%
% Material properties [infinitesimal deformation theory]
% MAT:0(Linear elastic) / MAT:1 (Elastoplastic) / MAT:2 (CMSG plasticity)
% PROP=[YOUNG NU] (MAT=0)
% PROP=[YOUNG NU SYIELD E0 N EDOT0 M] (MAT=1)
% PROP=[YOUNG NU H SYIELD] (MAT=2) 
% PROP=[YOUNG NU SYIELD ELLE N] (MAT=3)
%MAT=0;
%PROP=[260000 0.3];
%MAT=1;
%PROP=[260000, 0.3, 200, 0.2 10 0.1 20];
%MAT=2;
%PROP=[260000 0.3 1000 200];
MAT=3;
PROP=[260000 0.3 200 0.0 0.2];
% Prescribed displacements [Node, DOF, Value]
if Q==8 %[QUADRATIC ELEMENTS]
Dnodes=[1 2 0;10 2 0;2 2 0;14 2 0;3 2 0;  % u2=0 in the bottom line
        3 1 0;   % u1=0 in the bottom-right node (rigid body motion)
        7 2 0.3;18 2 0.3;8 2 0.3;21 2 0.3;9 2 0.3]; % u2=0.3 in the top line
else % Q=4 [LINEAR ELEMENTS]
Dnodes=[1 2 0;2 2 0;3 2 0;  % u2=0 in the bottom line
        3 1 0;   % u1=0 in the bottom-right node (rigid body motion)
        7 2 0.3;8 2 0.3;9 2 0.3]; % u2=0.3 in the top line  
end
%
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[10, 1E-6, 30];    
% Calling main function 
EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
