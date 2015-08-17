  function plotmesh(coords,connect,n,nelem,color)
% Function to plot the initial mesh.  
   f2D_4 = [1,2,3,4];
   f2D_8 = [1,5,2,6,3,7,4,8];
   hold on

       for lmn = 1:nelem
           for i = 1:n
               x(i,1:2) = coords(1:2,connect(i,lmn));
           end
           scatter(x(:,1),x(:,2),'MarkerFaceColor','r');
           if (n==4)
               patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor',color);
           elseif (n==8 || n==9)
               patch('Vertices',x,'Faces',f2D_8,'FaceColor','none','EdgeColor',color);
           end
       end
   axis equal
   hold off
end