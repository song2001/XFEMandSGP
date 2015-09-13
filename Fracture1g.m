% Fracture problem 1: plate loaded in tension
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
close all
clear variables
% Definition of global variables
global numcrack bcNodes edgNodes elemType stress_pnt
%
% Introduce manually the order of the interpolation (FEM X-FEM):
Q=[4 4];
if Q(1)==4
elemType = 'Q4' ;
elseif Q(1)==8
elemType = 'Q8' ;
end
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,XY1,LE,LE1]=ABAQUSmesh2D(Q);
% Identify boundary nodes
[bcNodes,edgNodes] = IdenBoundNodes(XY,LE) ;
% Crack definition
xCr(1).coor = [-14.1 0.0;0.0 0.0];
numcrack = size(xCr,2) ;
% enrichment type...enrType = geom (geometrical/fixed area enrichment), 
% enrType = ngeom (topological enrichment)
enrType = 'ngeom' ;
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
%MAT=0;
%PROP=[260000, 0.3];
%MAT=2;
%PROP=[260000 0.3 1000 200];
 MAT=3;
 PROP=[260000 0.3 200 0 0.2];
%
% Prescribed displacements [Node, DOF, Value]
% Introduce the vertical displacement of the upper edge
KLOAD=0.0011;
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
Dnodes(j,:)=[1 1 0];
%  
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[20, 1E-9, 50];
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA,stress_val]=EMPXFEM(xCr,XY,XY1,LE,LE1,MAT,PROP,Dnodes',SOL,Dloads',enrType,Q);

% Map results (Delaunay triangles and linear interpolation)
coord=[0.0001 0; 0.0002 0; 0.0003 0; 0.0004 0; 0.0005 0; 0.0006 0; 0.0007 0; 0.0008 0; 0.0009 0; 0.001 0;
       0.002 0; 0.003 0; 0.004 0; 0.005 0; 0.006 0; 0.007 0; 0.008 0; 0.009 0; 0.01 0; 0.015 0; 0.02 0; 0.025 0;
       0.03 0; 0.035 0; 0.04 0; 0.05 0; 0.06 0; 0.07 0; 0.08 0; 0.09 0; 0.1 0; 0.12 0; 0.14 0; 0.16 0; 0.18 0;
       0.2 0; 0.24 0; 0.28 0; 0.32 0; 0.38 0; 0.44 0; 0.5 0; 0.6 0; 0.7 0; 0.8 0; 0.9 0; 1 0; 1.5 0; 2 0; 3 0;
       5 0; 10 0; 19 0];
vq = griddata(stress_pnt(:,1),stress_pnt(:,2),stress_val,coord(:,1),coord(:,2), 'linear');
