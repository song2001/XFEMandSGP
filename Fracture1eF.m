% Fracture problem 1: plate loaded in tension
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
close all
clear variables
% Definition of global variables
global numcrack bcNodes edgNodes elemType
%
% Introduce manually the number of nodes per element:
Q=4;
if Q==4
elemType = 'Q4' ;
elseif Q==8
elemType = 'Q8' ;
end
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=ABAQUSmesh2D(Q);
% Identify boundary nodes
[bcNodes,edgNodes] = IdenBoundNodes(XY,LE) ;
% Crack definition
xCr(1).coor = [-14.1 0.0;0.0 0.0];
numcrack = size(xCr,2) ;
%plot the mesh before proceeding
plotNode = 'no' ;
plotMesh2(XY,LE,'b-',plotNode)
%crack plot
for k=1:size(xCr,2)
        for kj = 1:size(xCr(k).coor,1)-1
            cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            set(cr,'LineWidth',3);
        end
        for kj = 1:size(xCr(k).coor,1)
            plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        end
end
%
% Material properties
% MAT:0(Linear elastic) / MAT:1 (Elastoplastic) / MAT:2 (CMSG plasticity)
% PROP=[YOUNG NU] (MAT=0)
% PROP=[YOUNG NU SYIELD] (MAT=1) 
% PROP=[YOUNG NU SYIELD ELLE N] (MAT=2)
MAT=0;
PROP=[260000, 0.3];
%MAT=2;
%PROP=[260000 0.3 1000 200];
% MAT=3;
% PROP=[260000 0.3 200 0.0 0.2];
%
% Prescribed displacements [Node, DOF, Value]
% Constraint rigid body motion
Dnodes=[495 1 0];
%  
% Applied force [Element, Face, Traction component 1, Traction component 2]
% Face convention for quadrilateral elements: bottom - 1, right edge - 2, top - 3, left - 4.
Dloads=[918 1 0 -1; 917 1 0 -1; 916 1 0 -1; 915 1 0 -1; 914 1 0 -1; 913 1 0 -1; 912 1 0 -1;
    911 1 0 -1; 910 1 0 -1; 909 1 0 -1; 908 1 0 -1; 907 1 0 -1; 906 1 0 -1; 905 1 0 -1; 
    904 1 0 -1; 903 1 0 -1; 902 1 0 -1; 901 1 0 -1; % Load bottom part
    18 1 0 1; 17 1 0 1; 16 1 0 1; 15 1 0 1; 14 1 0 1; 13 1 0 1; 12 1 0 1;
    11 1 0 1; 10 1 0 1; 9 1 0 1; 8 1 0 1; 7 1 0 1; 6 1 0 1; 5 1 0 1; 
    4 1 0 1; 3 1 0 1; 2 1 0 1; 1 1 0 1]; % Load top part
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[1, 1E-9, 50];
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA]=EMPXFEM(xCr,XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
% Write to file
outfile2=fopen('S22.txt','w');
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,461));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,461));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,460));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,460));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,459));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,459));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,458));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,458));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,457));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,457));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,456));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,456));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,455));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,455));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,454));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,454));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,453));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,453));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,452));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,452));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,451));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,451));
fclose(outfile2);