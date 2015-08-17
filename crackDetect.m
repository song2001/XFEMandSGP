function [enrdomain]=crackDetect(Cr,enrdomain,node,element)

%Purpose
%identify elements that are in close proximity to the crack tip
 for icrack=1:size(Cr,2)
  Tip1 = Cr(icrack).coor(1,:) ;
  Tip2 = Cr(icrack).coor(size(Cr(icrack).coor,1),:) ;
  [enrdomain1,~] = ENRdomainf(Tip1,Tip2,node,element) ;
  [enrdomain2,~] = ENRdomainf(Tip2,Tip1,node,element) ;
  enrdomain = [enrdomain; enrdomain1; enrdomain2] ;
 end
enrdomain = unique(enrdomain) ;
end

