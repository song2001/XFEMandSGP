%
%================= ELEMENT STIFFNESS MATRIX ================================
%
function dsde = ElementK(dt,nelnodes,coord,materialprops,stress0,eplas0,MAT,rE,egrad0,ee0el,nelem,intpt,strain,N)

%
%  Assemble the element stiffness
%
%    Arguments;
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures
%      displacement(i,a)  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness
%
%  Set up integration points && weights    
%
   stressi0 = zeros(4,1);
   ee0=zeros(4,1);
      
      strain(4)=0; % Plane strain
%
        for j = 1 : 4
          stressi0(j) = stress0(j,intpt);
          ee0(j)= ee0el(j,intpt);
        end
      
%     Compute the material tangent stiffness (d stress/d strain)
%
   if MAT == 2
       for j=1:3
        deps(j) = strain(j);
       end  
    ep0 = eplas0(intpt);
    r=rE(intpt);
    D = ELAST(materialprops);
    dsde=combHardTanS(materialprops,D,deps',stressi0,r);
   elseif MAT==3
      gausscord(N,coord,nelem,intpt,nelnodes); 
    for j=1:3
      deps(j) = strain(j);
    end  
    ep0 = eplas0(intpt);
    eta0=egrad0(intpt);
    D = ELAST(materialprops);
    dsde=combHardTanC(materialprops,D,deps',ep0,ee0,eta0);
   elseif MAT==1
     ep0 = eplas0(intpt);
     dep = deplas(dt,stressi0,ep0,strain,materialprops);
     stress = materialstressP(stressi0,dep,strain,materialprops);
     dsde = materialstiffness(stress,ep0,dep,strain,materialprops);
   else
    dsde = ELAST(materialprops);
   end
end