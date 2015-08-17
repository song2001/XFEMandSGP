% Fracture problem 1: plate loaded in tension
%
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es
%
% Introduce manually the number of nodes per element:
Q=8;
% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=ABAQUSmesh2D(Q);
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
% PROP=[260000 0.3 200 0.00353 0.2];
%
% Prescribed displacements [Node, DOF, Value]
% Introduce the vertical displacement of the upper edge
KLOAD=0.0011;
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

    if XY(i,2)>49.999
    Dnodes(j,:)=[i 2 KLOAD];
    j=j+1;
    end
  end
  
% Constraint rigid body motion

    Dnodes(j,:)=[3 1 0];
%  
% Applied force [Element, Face, Traction component 1, Traction component 2]
Dloads=[];
% Solution controls [Number of steps, Convergence tolerance, maximum number
% of iterations]
SOL=[1, 1E-9, 50];    
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA]=EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
% Write to file
outfile2=fopen('S22.txt','w');
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6263));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6263));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6243));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6243));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6223));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6223));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6203));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6203));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6183));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6183));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6163));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6163));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6143));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6143));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6123));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6123));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6103));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6103));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6083));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6083));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6063));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6063));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6043));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6043));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6023));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6023));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6003));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6003));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5983));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5983));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5963));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5963));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5943));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5943));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5923));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5923));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5903));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5903));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5883));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5883));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5863));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5863));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5843));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5843));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5823));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5823));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5803));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5803));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5783));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5783));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5763));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5763));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5743));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5743));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5723));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5723));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5703));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5703));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5683));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5683));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,5283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,5283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2983));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2983));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2963));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2963));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2943));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2943));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2923));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2923));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2903));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2903));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2883));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2883));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2863));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2863));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2843));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2843));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2823));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2823));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2803));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2803));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2783));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2783));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2763));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2763));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2743));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2743));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2723));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2723));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2703));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2703));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2683));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2683));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2663));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2643));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2623));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2603));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2583));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2563));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2543));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2523));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2503));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2483));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2463));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2443));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2423));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2403));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2383));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2363));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2343));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2323));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2303));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2283));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2263));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2263));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2243));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2243));
fclose(outfile2);