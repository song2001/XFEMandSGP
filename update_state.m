function [stress_new,eplas_new,r_new,eta_new,ee_new,ep_new] = update_state(dt,ndof,coords,nelem,connect,materialprops,stress,eplas,dofs,MAT,rG,eta,ee,ep,nne,...
     enrich_node,elem_crk,type_elem,xTip,xVertex,split_elem,tip_elem,vertex_elem,pos,xCrk)

node=coords';
element=connect';
%   Assemble the global stiffness matrix
%
   lmncoord = zeros(ndof,nne);
   lmndof = zeros(ndof,nne);
%
%   Loop over all the elements
%
   for lmn = 1 : nelem
       
    sctr = element(lmn,:) ;
    
    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(lmn,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk,node,element) ;
        
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
      
      [lmnstress,lmneplas,lmnR,lmnETAP,lmnEE,lmnEPL] = update_el_state(dt,ndof,lmncoord,materialprops,lmnstress,lmneplas,U,MAT,lmn,lmnR,lmnETAP,lmnEE,lmnEPL,n,...
          type_elem,enrich_node,elem_crk,xVertex,node,element,W,Q,sctr,xCrk,tip_elem);
%
      if MAT > 0   
       for a = 1 : nintp
        for i = 1 : 3
          for j = 1 : 3
            stress_new(i,j,a,lmn) = lmnstress(i,j,a);
          end
        end
        for j = 1 : 4
            ee_new(j,a,lmn) = lmnEE(j,a);
            ep_new(j,a,lmn) = lmnEPL(j,a);
        end
        eplas_new(a,lmn) = lmneplas(a);
        r_new(a,lmn)=lmnR(a);
        eta_new(a,lmn)=lmnETAP(a);
       end    
      else
       for a = 1 : nintp
        for i = 1 : 2
          for j = 1 : 2
            stress_new(i,j,a,lmn) = stress(i,j,a,lmn)+lmnstress(i,j,a);
          end
        end
        for j = 1 : 4
            ee_new(j,a,lmn) = lmnEE(j,a);
            ep_new(j,a,lmn) = lmnEPL(j,a);
        end        
        eplas_new(a,lmn) = lmneplas(a);
        r_new(a,lmn)=lmnR(a);
        eta_new(a,lmn)=lmnETAP(a);
       end    
      end

  end
 
end