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
MAT=0;
PROP=[260000, 0.3];
%MAT=2;
%PROP=[260000 0.3 1000 200];
% MAT=3;
% PROP=[260000 0.3 200 0 0.2];
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
SOL=[1, 1E-9, 50];
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA,stress_pnt,stress_val]=EMPXFEM(xCr,XY,XY1,LE,LE1,MAT,PROP,Dnodes',SOL,Dloads',enrType,Q);

% Get results close to the crack tip
tri = delaunay(stress_pnt(:,1),stress_pnt(:,2));
coord=[0.0001 0; 0.0002 0; 0.0003 0; 0.0004 0; 0.0005 0; 0.0006 0; 0.0007 0; 0.0008 0; 0.0009 0; 0.001 0;
       0.002 0; 0.003 0; 0.004 0; 0.005 0; 0.006 0; 0.007 0; 0.008 0; 0.009 0; 0.01 0; 0.015 0; 0.02 0; 0.025 0;
       0.03 0; 0.035 0; 0.04 0; 0.05 0; 0.06 0; 0.07 0; 0.08 0; 0.09 0; 0.1 0; 0.12 0; 0.14 0; 0.16 0; 0.18 0;
       0.2 0; 0.24 0; 0.28 0; 0.32 0; 0.38 0; 0.44 0; 0.5 0; 0.6 0; 0.7 0; 0.8 0; 0.9 0; 1 0; 1.5 0; 2 0; 3 0;
       5 0; 10 0; 19 0];
   
%vq = griddata(x,y,v,xq,yq);
vq = griddata(stress_pnt(:,1),stress_pnt(:,2),stress_val,coord(:,1),coord(:,2), 'linear'); 
%In=0;   
%for i=1:size(coord(:,1))
%[N,dNdxi] = lagrange_basis(T3,coord(i,:));
% for j=1:size(tri(:,1))
%   for k=1:3
%   xv(k)=stress_pnt(tri(j,k),1);
%   yv(k)=stress_pnt(tri(j,k),2);
%   end  
% In=inpolygon(coord(i,1),coord(i,2),xv,yv);
%   if In==1
%   break
%   end
% end
%OBTAIN STRESSES XV AND XY CONTAIN THE "nODAL" COORDINATES WELL, WE ARE INTERESTED IN 
% OKAY, SOMEWHERE WE NEED TO INTRODUCE THE NATURAL COORDINATES, NOT THE
% ISOPARAMETRIC ONES
%end
% Write to file
%outfile2=fopen('S22.txt','w');
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,17,3159));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,5,3159));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,10,3159));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,9,3159));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3159));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,10,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,12,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,21,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,23,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,31,3041));
%fprintf(outfile2,'\n %9.10f',SIGMA(2,2,27,3041));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3042));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3042));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3043));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3043));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3044));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3044));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3045));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3045));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3046));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3046));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3047));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3047));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3048));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3048));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3049));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3049));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3050));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3050));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3051));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3051));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3052));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3052));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3053));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3053));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3054));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3054));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3055));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3055));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3056));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3056));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3057));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3057));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3058));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3058));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3059));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3059));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,4,3060));
% fprintf(outfile2,'\n %9.10f',SIGMA(2,2,2,3060));
% fclose(outfile2);

