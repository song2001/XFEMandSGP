% Edge crack under tension
% For validation purposes
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
close all
clear variables
% Definition of global variables
global xCr numcrack bcNodes edgNodes elemType
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
% Introduce the vertical displacement of the upper edge
KLOAD=0.005;
[NUMNP, ~] = size(XY);
j=1;
  for i=1:NUMNP
% 
% Impose load
%  
   if XY(i,2)<-49.999
    Dnodes(j,:)=[i 2 -KLOAD];
    j=j+1;
   end

   if XY(i,2)>49.999
    Dnodes(j,:)=[i 2 KLOAD];
    j=j+1;
   end
  end
  
% Constraint rigid body motion
Dnodes(j,:)=[495 1 0];
%  
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[1, 1E-9, 50];
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA]=EMPXFEM(xCr,XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
% Write to file
[Knum,Theta,xCrk] = KcalJint(xCrk,typeElem,Enrdomain,elemCrk,enrichNode,xVertex,vertexElem,pos,u,ipas,delta_inc,Knum,Theta,tipElem,splitElem) ;

a = 0.3 ;
C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
KAnalytical = C*P*sqrt(pi*a) 

Knumerical/KAnalytical