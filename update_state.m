function [stress_new,eplas_new,r_new,eta_new,ee_new,ep_new,stress_val] = update_state(dt,ndof,coords,nelem,connect,materialprops,stress,eplas,dofs,MAT,rG,eta,ee,ep,nne,...
     enrich_node,elem_crk,type_elem,xTip,xVertex,split_elem,tip_elem,vertex_elem,pos,xCrk,SOL,step)

 global STRAINP stress_pnt strainp1_val strainp2_val strainp3_val strainp4_val
node=coords';
element=connect';
%   Assemble the global stiffness matrix
%
   lmncoord = zeros(ndof,nne);
   lmndof = zeros(ndof,nne);
   
   stress_val = [ ] ;
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
%     STRAINP(lmn,:,1) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp1_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,2) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp2_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,3) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp3_val0,coordN(:,1),coordN(:,2), 'linear');
%     STRAINP(lmn,:,4) = griddata(stress_pnt(:,1),stress_pnt(:,2),strainp4_val0,coordN(:,1),coordN(:,2), 'linear');
    F1=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp1_val0);
    STRAINP(lmn,:,1) = F1(coordN);
    F2=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp2_val0);
    STRAINP(lmn,:,2) = F2(coordN); 
    F3=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp3_val0);
    STRAINP(lmn,:,2) = F3(coordN);  
    F4=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),strainp4_val0);
    STRAINP(lmn,:,2) = F4(coordN);        
    end
%Consider replacing griddata with scatteredinterpolant
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
      
      [lmnstress,lmneplas,lmnR,lmnETAP,lmnEE,lmnEPL,stress_val] = update_el_state(dt,ndof,lmncoord,materialprops,lmnstress,lmneplas,U,MAT,lmn,lmnR,lmnETAP,lmnEE,lmnEPL,n,...
          type_elem,enrich_node,elem_crk,xVertex,node,element,W,Q,sctr,xCrk,tip_elem,SOL,step,stress_val,coordN);
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

   if step==SOL(1) %Plot stress contours
      tri = delaunay(stress_pnt(:,1),stress_pnt(:,2)) ;    
      figure
      hold on
      for i=1:size(tri,1)
       fill(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:))) ;
      end
   title('XFEM Stress')
   shading interp
   colorbar 
   end

end