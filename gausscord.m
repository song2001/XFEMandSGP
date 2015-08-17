function gausscord(N,coord,nelem,intpt,n)
% Computation of the physical coordinates of the Gauss points
global CORDE

if n==8
cordX=N(1)*coord(1,1)+N(2)*coord(1,2)+N(3)*coord(1,3)+N(4)*coord(1,4)+N(5)*coord(1,5)+N(6)*coord(1,6)+N(7)*coord(1,7)+N(8)*coord(1,8);
cordY=N(1)*coord(2,1)+N(2)*coord(2,2)+N(3)*coord(2,3)+N(4)*coord(2,4)+N(5)*coord(2,5)+N(6)*coord(2,6)+N(7)*coord(2,7)+N(8)*coord(2,8);
else
cordX=N(1)*coord(1,1)+N(2)*coord(1,2)+N(3)*coord(1,3)+N(4)*coord(1,4);
cordY=N(1)*coord(2,1)+N(2)*coord(2,2)+N(3)*coord(2,3)+N(4)*coord(2,4);
end
CORDE(nelem,1,intpt)=cordX;
CORDE(nelem,2,intpt)=cordY;

end