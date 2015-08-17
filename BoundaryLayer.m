% Boundary layer formulation for the sharp crack problem
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
format long
% Introduce manually the number of nodes per element:
Q=4;
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=ABAQUSmesh2D(Q);
%
% Material properties
% MAT:0(Linear elastic) / MAT:1 (Elastoplastic) / MAT:2 (CMSG plasticity)
% PROP=[YOUNG NU] (MAT=0)
% PROP=[YOUNG NU SYIELD] (MAT=1) 
% PROP=[YOUNG NU SYIELD ELLE N] (MAT=2)
%MAT=0;
%PROP=[260000, 0.3];
MAT=2;
PROP=[260000 0.3 1000 200];
% MAT=3;
% PROP=[260000 0.3 200 0.00353 0.2];
%
% Prescribed displacements [Node, DOF, Value]
% Assign manually a value to K
KLOAD=205.57;
[NUMNP, ~] = size(XY);
j=1;
  for i=1:NUMNP
%
% Constraint vertical displacements
%

    if(XY(i,2)<0.000000001) && (XY(i,1)>=0)
    Dnodes(j,:)=[i 2 0];
    j=j+1;
    end 
% 
% Impose load
%    
   radi=sqrt(XY(i,2)*XY(i,2)+XY(i,1)*XY(i,1));
 
    if(radi>41)
        if(XY(i,1)>0)
        theta=atan(XY(i,2)/XY(i,1));
        else
        theta=pi-atan(XY(i,2)/abs(XY(i,1)));     
        end
    Dnodes(j,:)=[i 2 (1+PROP(2))*(KLOAD/PROP(1))*sqrt(radi/(2*pi))*(3-4*PROP(2)-cos(theta))*sin(theta/2)];
    Dnodes(j+1,:)=[i 1 (1+PROP(2))*(KLOAD/PROP(1))*sqrt(radi/(2*pi))*(3-4*PROP(2)-cos(theta))*cos(theta/2)];
    j=j+2;
    end
  end
%  
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[100, 1E-9, 50];    
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA]=EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
%
% We now obtain the opening stress distribution ahead of the crack-tip
for m=1:NUMNP
 XNELM(1,m)=0.0;
 XNELM(2,m)=0.0;
 SIGMANOD22(m)=0.0;
end
% Equivalence of (node) and Gauss point
NG(1)=1;
NG(2)=2;
NG(3)=4;
NG(4)=3;
[NE, ~] = size(LE);
for n=1:NE
 for u=1:4
 INOD=LE(n,u);
 XNELM(1,INOD)=XNELM(1,INOD)+SIGMA(2,2,NG(u),n);
 XNELM(2,INOD)=XNELM(2,INOD)+1;
 end 
end
for m=1:NUMNP
 SIGMANOD22(m)=XNELM(1,m)/XNELM(2,m);
end

% Write to file
outfile2=fopen('S22.txt','w');
for m=1:NUMNP
 if (XY(m,1)>0 && XY(m,2)<0.00000001)
 fprintf(outfile2,'\n %9.10f %9.10f %9.10f',XY(m,1),XY(m,2),SIGMANOD22(m));
 end   
end
fclose(outfile2); 

