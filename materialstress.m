%================= Material Stress ==================================
%
%   Computes stress sigma_{ij} given strain epsilon_{ij}
%
function stress = materialstress(ndof,strain,materialprops)

   C = ELAST(materialprops);
   stress = zeros(ndof,ndof);
   stressT=C*strain;
  stress(1,1) = stressT(1);
  stress(2,2) = stressT(2);
  stress(1,2) = stressT(3);
  stress(2,1) = stress(1,2);
  end