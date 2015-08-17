% Example of input file for the code EMPFEM
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
MAT=0;
PROP=[260000, 0.3];
%MAT=1;
%PROP=[260000, 0.3, 200, 0.2 10 0.1 20];
%MAT=2;
%PROP=[260000 0.3 1000 200];
%MAT=3;
%PROP=[260000 0.3 200 0 0.2];
%
% Prescribed displacements [Node, DOF, Value]
if Q==8 %[QUADRATIC ELEMENTS]
Dnodes=[1 1 0;21 1 0;2 1 0;25 1 0;3 1 0;28 1 0;4 1 0;31 1 0;5 1 0; % u1=0 in the bottom line
        1 2 0;21 2 0;2 2 0;25 2 0;3 2 0;28 2 0;4 2 0;31 2 0;5 2 0;   % u2=0 in the bottom line
        16 2 0;44 2 0;17 2 0;47 2 0;18 2 0;49 2 0;19 2 0;51 2 0;20 2 0;   % u2=0 in the top line
        16 1 1;44 1 1;17 1 1;47 1 1;18 1 1;49 1 1;19 1 1;51 1 1;20 1 1]; % u1=1 in the top line
else % Q=4 [LINEAR ELEMENTS]
Dnodes=[1 1 0;2 1 0;3 1 0;4 1 0;5 1 0; % u1=0 in the bottom line
        1 2 0;2 2 0;3 2 0;4 2 0;5 2 0;   % u2=0 in the bottom line
        16 2 0;17 2 0;18 2 0;19 2 0;20 2 0;   % u2=0 in the top line
        16 1 1;17 1 1;18 1 1;19 1 1;20 1 1]; % u1=1 in the top line    
end
%
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[10, 1E-3, 50];    
% Calling main function 
EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
