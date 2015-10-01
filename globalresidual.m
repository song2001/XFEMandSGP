function resid = globalresidual(dt,ndof,coords,nelem,connect,materialprops,stress,eplas,dofs,MAT,rG,ee,ep,eta,nne,tu,...
    enrich_node,elem_crk,type_elem,xTip,xVertex,split_elem,tip_elem,vertex_elem,pos,xCrk,step,nit)

global elemType STRAINP stress_pnt strainp1_val strainp2_val strainp3_val strainp4_val
node=coords';
element=connect';
%
%   Assemble the global stiffness matrix
%

   resid = zeros(tu,1);
   lmncoord = zeros(ndof,nne);
   lmndof = zeros(ndof,nne);
   rel = zeros(ndof*nne,ndof*nne);
   
   % Empty strainp1_val, make strainp1_valTemp=strainp1
   strainp1_val0=strainp1_val; strainp2_val0=strainp2_val; strainp3_val0=strainp3_val; strainp4_val0=strainp4_val;
   strainp1_val=[]; strainp2_val=[]; strainp3_val=[]; strainp4_val=[];
   
%
%   Loop over all the elements
%
   for lmn = 1 : nelem

    sctr = element(lmn,:) ; 
    coordN=node(sctr',:);
    [Nlines, ~] = size(coordN);
    if Nlines>4
     coordN(5,:)=[];  coordN(5,:)=[];  coordN(5,:)=[]; coordN(5,:)=[];
    end
    
    if MAT==3
    if step>1 || nit>2
%     STRAINP(lmn,:,1) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp1_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,2) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp2_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,3) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp3_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,4) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp4_val0,coordN(:,1),coordN(:,2), 'linear');
% vq = griddata(x,y,v,xq,yq) fits a surface of the form v = f(x,y) to the scattered data in the vectors (x,y,v). 
% The griddata function interpolates the surface at the query points specified by (xq,yq) and returns the interpolated values, 
% vq. The surface always passes through the data points defined by x and y.
    F1=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp1_val0);
    STRAINP(lmn,:,1) = F1(coordN);
    F2=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp2_val0);
    STRAINP(lmn,:,2) = F2(coordN); 
    F3=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp3_val0);
    STRAINP(lmn,:,2) = F3(coordN);  
    F4=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp4_val0);
    STRAINP(lmn,:,2) = F4(coordN);    
    
    end %Consider replacing griddata with scatteredinterpolant
    end
    
    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(lmn,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk,node,element);
        
% Kind of strain-displacement matrix to be computed / Global DOFs associated with the element  
    sctrB = [ ] ;
    for k = 1:size(xCrk,2)
        sctrB = [sctrB assembly(lmn,enrich_node(:,k),pos(:,k),k,element)] ;
    end    
    
% Computation of the displacement    
    U = [ ];
    for k = 1:size(xCrk,2)
        U = [U; element_disp(lmn,pos(:,k),enrich_node(:,k),dofs,k,element)];
    end
    
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1 : nne
        for i = 1 : ndof
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1 : ndof
          lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
      n = nne;
      nintp = size(W,1);
      lmnstress = zeros(4,nintp);
      lmneplas = zeros(nintp,1);
      lmnR=zeros(nintp,1);
      lmnEPL=zeros(4,nintp);
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
        lmnEPL(1,a)=ep(1,a,lmn);
        lmnEPL(2,a)=ep(2,a,lmn);
        lmnEPL(3,a)=ep(3,a,lmn);
        lmnEPL(4,a)=ep(4,a,lmn);
        lmnEE(1,a)=ee(1,a,lmn);
        lmnEE(2,a)=ee(2,a,lmn);
        lmnEE(3,a)=ee(3,a,lmn);
        lmnEE(4,a)=ee(4,a,lmn);
        lmnETAP(a)=eta(a,lmn);       
      end
      end
      
   if MAT > 0   
    strain = zeros(4,1);
   else
   strain = zeros(3,1);
   end

    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,lmn,type_elem,enrich_node(:,k),elem_crk,xVertex,k,node,element,MAT,tip_elem)] ;
        end
        Ppoint =  N' * node(sctr,:);
        strain = B*U ;   
            
      stressE = elresid(dt,ndof,lmncoord,materialprops,lmnstress,lmneplas,strain,MAT,lmn,lmnR,lmnETAP,lmnEE,lmnEPL,kk,coordN,Gpt);

      stresst(1)=stressE(1,1);
      stresst(2)=stressE(2,2);
      stresst(3)=stressE(1,2);
      elres=B'*stresst'*W(kk)*det(JO) ;
      
      resid(sctrB) = resid(sctrB) + elres ;
      
   end
 end
end