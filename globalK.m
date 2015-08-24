function GK=globalK(dt,ndof,coords,nelem,connect,materialprops,stress,eplas,w,MAT,rG,ee,eta,nne,tu,...
    enrich_node,elem_crk,type_elem,xTip,xVertex,split_elem,tip_elem,vertex_elem,pos,xCrk,step,nit)

%declare global variables here
global elemType stress_pnt

node=coords';
element=connect';

%
%   Assemble the global stiffness matrix
%
 
   GK = zeros(tu,tu);
   lmncoord = zeros(2,nne);
   lmndof = zeros(ndof,nne);
 
   q = [] ; % X-FEM feature
%
%   Loop over all the elements
%
   for lmn = 1:nelem
       
   sctr = element(lmn,:);  
   
%choose Gauss quadrature rules for elements
   [W,Q]=gauss_rule(lmn,enrich_node,elem_crk,xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk,node,element);   
 
%Transfrom these Gauss points to global coords for plotting ONLY
    for igp = 1:size(W,1)
        gpnt = Q(igp,:) ;
        [N,~] = lagrange_basis(elemType,gpnt) ;
        pt = N' * node(sctr,:);
        
        if step==1 && nit==1
          stress_pnt = [stress_pnt; pt] ;  
        end
        
%         if lmn==3041
%         GpTT(igp,:)=Gpnt;    
%         end
%       node1=[0.5774 0.5774; 0.5774 -0.5774;-0.5774 0.5774;-0.5774 -0.5774];
%        %node1=[1.4226 1.4226; -1.4226 1.4226;-1.4226 -1.4226;1.4226 -1.4226]; % For lmn==3159
%        node1=[-1.4226 -1.4226; 1.4226 -1.4226;1.4226 1.4226;-1.4226 1.4226]; % For lmn==3041
%        if lmn==3159 || lmn==3041
%        Temp= N'*node1;    
%        INTP(lmn,igp,1)= Temp(1);
%        INTP(lmn,igp,2)= Temp(2);
%        if lmn==3159
%        Prueba(igp,:)=[Temp];
%        end
%        if lmn==3041
%        Prueba1(igp,:)=[Temp];
%        end       
%        end
%        if lmn==3041    
%        Temp= N'*node1;
%        INTP(lmn,igp,1)= Temp(1);
%        INTP(lmn,igp,2)= Temp(2);
%        end       
        q = [q;pt] ;
    end   
    
% Kind of strain-displacement matrix to be computed / Global DOFs associated with the element  
    sctrB = [ ] ;
    for k = 1:size(xCrk,2)
        sctrB = [sctrB assembly(lmn,enrich_node(:,k),pos(:,k),k,element)] ;
    end    
    
% Computation of the displacement    
    U = [ ];
    for k = 1:size(xCrk,2)
        U = [U; element_disp(lmn,pos(:,k),enrich_node(:,k),w,k,element)];
    end   
    
%%% THIS COMMES FROM EMPFEM (necessary to compute the constitutive matrix)
     for a = 1:nne
       for i = 1:2
         lmncoord(i,a) = coords(i,connect(a,lmn));
       end
       for i = 1:ndof
         lmndof(i,a) = w(ndof*(connect(a,lmn)-1)+i);
       end
     end
     n = nne;
     nintp = size(W,1);
     lmnstress = zeros(4,nintp);
     lmneplas = zeros(nintp,1);
     lmnR=zeros(nintp,1);
     lmnEE=zeros(4,nintp);
     lmnETAP=zeros(nintp,1);
     if MAT > 0
      for a = 1 : nintp
        lmnstress(1,a) = stress(1,1,a,lmn);
        lmnstress(2,a) = stress(2,2,a,lmn);
        lmnstress(4,a) = stress(3,3,a,lmn);
        lmnstress(3,a) = stress(1,2,a,lmn);
        lmneplas(a) = eplas(a,lmn);
        lmnR(a)=rG(a,lmn);
        lmnEE(1,a)=ee(1,a,lmn);
        lmnEE(2,a)=ee(2,a,lmn);
        lmnEE(3,a)=ee(3,a,lmn);
        lmnEE(4,a)=ee(4,a,lmn);
        lmnETAP(a)=eta(a,lmn);
      end  
     end    

%%% THIS COMES FROM EMPFEM     
     
    %loop over Gauss points
    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,lmn,type_elem,enrich_node(:,k),elem_crk,xVertex,k,node,element,MAT,tip_elem)] ;
        end
        Ppoint =  N' * node(sctr,:);
        eps_sub = B*U ;
        
        % We need to compute C here and for that we need to call the FEM
        % part
        
        C = ElementK(dt,n,lmncoord,materialprops,lmnstress,lmneplas,MAT,lmnR,lmnETAP,lmnEE,lmn,kk,eps_sub,N);
            
        GK(sctrB,sctrB) = GK(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
    end %Gauss points   
   end %Elements
end %Function