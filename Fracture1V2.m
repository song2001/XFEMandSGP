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
%MAT=0;
%PROP=[260000, 0.3];
%MAT=2;
%PROP=[260000 0.3 1000 200];
 MAT=3;
 PROP=[260000 0.3 200 0.005 0.2];
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
SOL=[100, 1E-9, 50];    
% Calling main function and storing the stress tensor for every integration
% point and every element
[SIGMA]=EMPFEM(XY,LE,MAT,PROP,Dnodes',SOL,Dloads');
% Write to file
outfile2=fopen('S22.txt','w');
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,8006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,8006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7386));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7386));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7366));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7366));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7346));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7346));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7326));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7326));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7306));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7306));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7286));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7286));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7266));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7266));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7246));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7246));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7226));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7226));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7206));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7206));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,7006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,7006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,6406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,6406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2586));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2566));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2546));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2526));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2506));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2486));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2466));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2446));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2426));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2406));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2386));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2386));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2366));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2366));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2346));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2346));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2326));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2326));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2306));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2306));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2286));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2286));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2266));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2266));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2246));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2246));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2226));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2226));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2206));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2206));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2186));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2166));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2146));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2126));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2106));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2086));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2066));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2046));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2026));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,2006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,2006));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1986));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1966));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1946));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1926));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1906));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1886));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1866));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1846));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1826));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1806));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1786));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1766));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1746));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1726));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1706));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1686));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1666));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1646));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1626));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,3,1606));
fprintf(outfile2,'\n %9.10f',SIGMA(2,2,1,1606));
fclose(outfile2);