function xi = integrationpoints(ncoord,nelnodes,npoints)

   xi = zeros(ncoord,npoints);
%
%  1D elements
%
   if (ncoord == 1) 
     if (npoints==1) 
       xi(1,1) = 0.;
     elseif (npoints == 2) 
       xi(1,1) = -0.5773502692;
       xi(1,2) = -xi(1,1);
     elseif (npoints == 3) 
       xi(1,1) = -0.7745966692;
       xi(1,2) = 0.0;
       xi(1,3) = -xi(1,1);
     end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 || nelnodes == 6 ) 
       if (npoints == 1) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
       elseif (npoints == 3) 
         xi(1,1) = 0.6;
         xi(2,1) = 0.2;
         xi(1,2) = 0.2;
         xi(2,2) = 0.6;
         xi(1,3) = 0.2;
         xi(2,3) = 0.2;
       elseif (npoints == 4) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
         xi(1,2) = 0.6;
         xi(2,2) = 0.2;
         xi(1,3) = 0.2;
         xi(2,3) = 0.6;
         xi(1,4) = 0.2;
         xi(2,4) = 0.2;
       end
%
%    Rectangular element
%                  
     elseif ( nelnodes==4 || nelnodes==8 ) 

       if (npoints == 1) 
         xi(1,1) = 0.;
         xi(2,1) = 0.;
       elseif (npoints == 4) 
         xi(1,1) = -0.5773502692;
         xi(2,1) = xi(1,1);
         xi(1,2) = -xi(1,1);
         xi(2,2) = xi(1,1);
         xi(1,3) = xi(1,1);
         xi(2,3) = -xi(1,1);
         xi(1,4) = -xi(1,1);
         xi(2,4) = -xi(1,1);
       elseif (npoints == 9) 
         xi(1,1) = -0.7745966692;
         xi(2,1) = xi(1,1);
         xi(1,2) = 0.0;
         xi(2,2) = xi(1,1);
         xi(1,3) = -xi(1,1);
         xi(2,3) = xi(1,1);
         xi(1,4) = xi(1,1);
         xi(2,4) = 0.0;
         xi(1,5) = 0.0;
         xi(2,5) = 0.0;
         xi(1,6) = -xi(1,1);
         xi(2,6) = 0.0;
         xi(1,7) = xi(1,1);
         xi(2,7) = -xi(1,1);
         xi(1,8) = 0.;
         xi(2,8) = -xi(1,1);
         xi(1,9) = -xi(1,1);
         xi(2,9) = -xi(1,1);
       end
    end
%
end