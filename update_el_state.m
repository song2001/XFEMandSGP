function [stress1,ep1,r1,etap1,ee1,epl1,stress_pnt,stress_val] = update_el_state(dt,ndof,coord,materialprops,stress0,eplas0,U,MAT,nelem,rE,egrad0,ee0el,epl0el,n,...
    type_elem,enrich_node,elem_crk,xVertex,node,element,W,Q,sctr,xCrk,tip_elem,SOL,step,stress_pnt,stress_val)

global elemType
%
%  Compute the stress and accumulated plastic strain at the end of a load
%  increment at all integration points in an element
%
%    Arguments:
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords[i,a]        ith coord of ath node
%      materialprops      Material properties passed on to constitutive procedures
%      displacement[i,a]  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi[i,inpt]         ith local coord of integration point no. intpt
%      w[intpt]           weight for integration point no. intpt
%      N[a]               Shape function associated with ath node on element
%      dNdxi[a,i]         Derivative of ath shape function wrt ith local coord
%      dNdx[a,i]          Derivative of ath shape function wrt ith global coord
%      dxdxi[i,j]         Derivative of ith global coord wrt jth local coord
%      dxidx[i,j]         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain[i,j]        strain_ij components
%      stress[i,j]        stress_ij components
%      r[row]             Residual vector
%
%
   npoints = 4;
   dxdxi = zeros(ndof,ndof);
   dxidx = zeros(ndof,ndof);
   if MAT > 0   
    strain = zeros(4,1);
   else
   strain = zeros(3,1);
   end
   stressi0 = zeros(4,1);
   ee0=zeros(4,1);
   epl0=zeros(4,1);
   
%
%  Set up integration points and weights    
%
   xilist = integrationpoints(ndof,n,npoints);
   w = integrationweights(ndof,n,npoints);   
   %
%  Loop over the integration points
%
%   for intpt = 1 : npoints
       
    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        pt = N' * node(sctr,:);
        Gpnt = N'*node(sctr,:) ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,nelem,type_elem,enrich_node(:,k),elem_crk,xVertex,k,node,element,MAT,tip_elem)] ;
        end
        Ppoint =  N' * node(sctr,:);
        strain = B*U ;       
    
      if MAT > 0
      strain(4)=0.0; % Plane strain
      ep0 = eplas0(kk);
      eta0=egrad0(kk);
      end      
%
%     Compute the stress
%
        for j = 1 : 4
          stressi0(j) = stress0(j,kk);
          ee0(j)= ee0el(j,kk);
          epl0(j)= epl0el(j,kk);
        end
%     Compute the material tangent stiffness (d stress/d strain)
%
      if MAT == 2 
       for j=1:3
        deps(j) = strain(j);
       end
       r=rE(kk);
       D = ELAST(materialprops);
       [stress, dep,dr]=combHardS(materialprops,D,deps',stressi0,r);   
      elseif MAT==3
       for j=1:3
        deps(j) = strain(j);
       end 
       [stress, dep, dee, depl, detap]=combHardC(materialprops,deps',ep0,ee0,epl0,eta0,kk,nelem);
      elseif MAT==1
        dep = deplas(dt,stressi0,ep0,strain,materialprops);
        stress = materialstressP(stressi0,dep,strain,materialprops);
      else
       stress = materialstress(ndof,strain,materialprops);
      end
      
      if step==SOL(1) %Plot stress contours
        stress_pnt = [stress_pnt; pt] ;
        stress_val = [stress_val; stress(2,2)] ;      
      end    
      
      if MAT > 0   
       for i = 1 : 3
         for j = 1 : 3
           stress1(i,j,kk) = stress(i,j);
         end
       end
      else
       ep0=0;
       dep=0;
       dr=0;
       dee=0;
       depl=0;
       detap=0;
       for i = 1 : 2
        for j = 1 : 2
          stress1(i,j,kk) = stress(i,j);
        end
       end
      end

     ep1(kk) = ep0 + dep;
     if MAT == 3
     dr=0;
     etap1(kk)=detap;
     ee1(:,kk)=dee;
     epl1(:,kk)=depl;
     else
     etap1=zeros(kk); ee1=zeros(4,kk); epl1=zeros(4,kk);
     end
     r1(kk)= dr;

      
   end

end