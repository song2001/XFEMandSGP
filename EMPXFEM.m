function [SIGMA,stress_val]=EMPXFEM(xCrk,XY,XY1,LE,LE1,MAT,PROP,Dnodes,SOL,Dloads,enrType,Q)

global STRAINP strainp1_val strainp2_val strainp3_val strainp4_val stress_pnt
%
% 2D plane strain X-FEM code
% Current implementation: linear and quadratic elements with 4 Gauss points
% Material models: elasticity, plasticity, CMSG plasticity
% Displacement controlled models (currently)
% Emilio Martínez-Pañeda
% martinezemilio@uniovi.es

[NNodes, Ndof]=size(XY);
[NE, NNE] = size(LE);
[~, nfix] = size(Dnodes);
[~, ndload] = size(Dloads);
%plotmesh(XY',LE',NNE,NE,'g');

Enrdomain = [ ] ; tipElem = [ ] ; splitElem = [ ] ; vertexElem = [ ] ;
% Identify enriched elements
if Q(1)==8 && Q(2)==4
[Enrdomain] = crackDetect(xCrk,Enrdomain,XY1,LE1) ;
% Classify elements
[typeElem,elemCrk,tipElem,splitElem,vertexElem,xTip,xVertex,enrichNode]=nodeDetect(xCrk,Enrdomain,XY1,LE1) ;
enrichNode(numel(XY(:,1))) = 0;
else
[Enrdomain] = crackDetect(xCrk,Enrdomain,XY,LE) ;
% Classify elements
[typeElem,elemCrk,tipElem,splitElem,vertexElem,xTip,xVertex,enrichNode]=nodeDetect(xCrk,Enrdomain,XY,LE) ;
end

if strcmp(enrType,'geom')
  % find the area of tip element...
  econ = LE(tipElem,:) ;
  nds = XY(econ,:) ;
  [geominfo,cpmo,iner] = polygeom(nds(:,1),nds(:,2)) ;
  tipArea = geominfo(1) ;
  % geometrix/fixed area enrichment
  cnt = 0 ;
  fa_fac = 50; % fixed area factor..!!
  Renr = fa_fac*tipArea ;
  %Renr = 1 ;
  xcent = xCrk(1).coor(2,:) ;
  figure(1)
  theta = -pi:0.1:pi ;
  xo = xcent(1) + Renr*cos(theta) ;
  yo = xcent(2) + Renr*sin(theta) ;
  plot(xo,yo,'k-')
  for in = 1:NNodes
    x = XY(in,:) ;
    dist = sqrt( (x(1)-xcent(1))^2+(x(2) - xcent(2))^2) - Renr ;
       if dist < 0
         cnt = cnt + 1 ;
         enrichNode(in,1) = 1 ;
       end
   end
end
% Initialize stiffness matrix and force vector
% Take into account that each split node is enriched by ONE function, H(x)
% while each tip node is enriched by FOUR functions, B_i(x), i=1,2,3,4
% total dof = numnode*ndof + numsplitnode*1*ndof + numtipnode*4*ndof

    split = 0 ; tip = 0 ;
    for k=1:size(xCrk,2)
        split = split + size(find(enrichNode(:,k) == 2), 1) ;
        tip = tip + size(find(enrichNode(:,k) == 1),1 ) ;
    end

    pos = posi(xCrk,NNodes,enrichNode) ;    
    
    TU = NNodes*Ndof + split*1*Ndof + tip*4*Ndof ;

% Plot enriched nodes    
for k=1:size(xCrk,2)
 for kj = 1:size(xCrk(k).coor,1)-1
  cr = plot(xCrk(k).coor(kj:kj+1,1),xCrk(k).coor(kj:kj+1,2),'r-') ;
  set(cr,'LineWidth',3);
 end
 for kj = 1:size(xCrk(k).coor,1)
   plot(xCrk(k).coor(kj,1),xCrk(k).coor(kj,2),'ro','MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
 end
   split_nodes = find(enrichNode(:,k) == 2);
   tip_nodes   = find(enrichNode(:,k) == 1);
   n1 = plot(XY(split_nodes,1),XY(split_nodes,2),'r*');
   n2 = plot(XY(tip_nodes,1),XY(tip_nodes,2),'rs');
   set(n1,'MarkerSize',15);
   set(n2,'MarkerSize',15);
end
    
% Initialize the displacement and the incremental displacements
w = zeros(TU,1);
u=w;

% Initialize history dependent variables

SIGMA = zeros(3,3,300,NE);
EP=zeros(300,NE);
rG=zeros(300,NE);
EE=zeros(4,300,NE);
EPL=zeros(4,300,NE);
ETAP=zeros(300,NE);
strainp1_val=[]; 
strainp2_val=[];
strainp3_val=[];
strainp4_val=[];
stress_pnt=[];

STRAINP=zeros(NE,4,4);

% Store solver controls

dt = 2./SOL(1); 
forcevdisp = zeros(2,SOL(1)+1);
forcevdisp(1,1) = 0;
forcevdisp(2,1) = 0;


for step = 1 : SOL(1)

 loadfactor = step/SOL(1);
 err1 = 1.;    % Error/residual
 nit = 0;      % Number of iterations
 
 fprintf(1,'\n Step %f Load %f\n',step,loadfactor);
 
 while ((err1>SOL(2)) && (nit<SOL(3)))          % Newton Raphson loop
   
   nit = nit + 1;
   % Compute the global stiffness matrix, force vector and residual   
   GK=globalK(dt,Ndof,XY',NE,LE',PROP,SIGMA,EP,w,MAT,rG,EE,ETAP,NNE,TU,enrichNode,elemCrk,...
      typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,pos,xCrk,step,nit); 
   GF=globaltraction(Ndof,ndload,XY',LE',Dloads,w,NNE,TU);
   R = globalresidual(dt,Ndof,XY',NE,LE',PROP,SIGMA,EP,w,MAT,rG,EE,EPL,ETAP,NNE,TU,...
       enrichNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,pos,xCrk,step);
   b = loadfactor*GF - R;
  
   % Prescribed displacements
         
         for n = 1 : nfix
           rw = Ndof*(Dnodes(1,n)-1) + Dnodes(2,n);
           for cl = 1 : TU %Should I have changed from Ndof*NNodes?
             GK(rw,cl) = 0;
           end
           GK(rw,rw) = 1.;
           b(rw) = loadfactor*Dnodes(3,n)-w(rw)-u(rw);
         end
 
   % Solve system

         dw = GK\b;

   % Check convergence

         w = w + dw;
         wnorm = dot(w,w);
         err1 = dot(dw,dw);
         err2 = dot(b,b);
         err1 = sqrt(err1/wnorm);
         err2 = sqrt(err2)/(Ndof*NNodes);
         fprintf(1,'Iteration number %d Correction %f Residual %f tolerance %f\n',nit,err1,err2,SOL(2));
         if nit == SOL(3)
         fprintf(1,'****** NEWTON-RAPHSON DID NOT CONVERGE');
         stop
         end    
  end
    
% Update the stress and accumulated plastic strain
     [SIGMA,EP,rG,ETAP,EE,EPL,stress_val] = update_state(dt,Ndof,XY',NE,LE',PROP,SIGMA,EP,w,MAT,rG,ETAP,EE,EPL,NNE,...
         enrichNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,pos,xCrk,SOL,step);
     
% Update the total displacecment
     u = u + w;
%
end

% 
% Print the solution
%
outfile=fopen('FEM_results.txt','w');
fprintf(outfile,'\n\n Step %f Load %f\n',step,loadfactor);

print_results(outfile,Ndof,NNodes,XY',NE,LE',SIGMA,EP,u,NNE);

fclose(outfile);
end 