% Fracture problem 1: plate loaded in tension
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
close all
clear variables
% Definition of global variables
global numcrack bcNodes edgNodes elemType %Equiv
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
[SIGMA]=EMPXFEM(xCr,XY,LE,MAT,PROP,Dnodes',SOL,Dloads',enrType);
% Write to file
outfile2=fopen('S22.txt','w');
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,17,3159));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,5,3159));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,10,3159));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,9,3159));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3159));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,10,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,12,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,21,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,23,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,31,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,27,3041));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3042));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3042));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3043));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3043));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3044));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3044));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3045));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3045));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3047));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3047));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3048));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3048));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3049));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3049));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3050));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3050));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3051));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3051));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3052));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3052));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3053));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3053));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3054));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3054));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3055));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3055));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3056));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3056));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3057));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3057));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3058));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3058));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3059));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3059));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3060));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3060));
fclose(outfile2);

