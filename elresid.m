function stress = elresid(dt,ndof,coord,materialprops,stress0,eplas0,strain,MAT,nelem,rE,egrad0,ee0el,epl0el,kk,coordN,Gpt)


%
%  Assemble the element residual force
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

      if MAT > 0
      strain(4)=0.0; % Plane strain
      ep0 = eplas0(kk);
      end
   stressi0 = zeros(4,1);
   ee0=zeros(4,1);
   epl0=zeros(4,1);
    

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
       [stress, ~]=combHardS(materialprops,D,deps',stressi0,r);
      elseif MAT==3
       for j=1:3
        deps(j) = strain(j);
       end
       eta0=egrad0(kk);
       [stress,~,~,~,~]=combHardC(materialprops,deps',ep0,ee0,epl0,eta0,nelem,coordN,Gpt(1),Gpt(2));
      elseif MAT==1
       dep = deplas(dt,stressi0,ep0,strain,materialprops);
       stress = materialstressP(stressi0,dep,strain,materialprops);    
      else
       stress = materialstress(ndof,strain,materialprops);
      end
  
end