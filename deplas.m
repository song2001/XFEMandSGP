%======================= deplas ===================
%
%  Computes plastic strain increment de_e given strain increment
%
function e = deplas(dt,stress,eplas,dstrain,materialprops)

%  Youngs modulus and Poissons ratio
   E = materialprops(1);
   nu = materialprops(2);
   Y = materialprops(3);
   e0 = materialprops(4);
   n = materialprops(5);
   edot0 = materialprops(6);
   m = materialprops(7);
%
%  S is the deviatoric stress predictor
%
   S = zeros(3,3);
   de = zeros(3,3);
   dl = [ [1,0,0];[0,1,0];[0,0,1] ];  

   devol = trace(dstrain);
   p = trace(stress);
   
   sequiv = 0.;
   for i = 1 : 3
     for j = 1 : 3
        de(i,j) = dstrain(i,j) - dl(i,j)*devol/3;
        S(i,j) = stress(i,j) - dl(i,j)*p/3 + E/(1+nu)*de(i,j);
        sequiv = sequiv + S(i,j)*S(i,j);
     end
   end
   sequiv = sqrt(1.5*sequiv);

   e = 10^(-15);
   err = Y;
   tol = 10^(-06)*Y;
   if (sequiv*edot0 == 0)
     e = 0.;
   else
     while (err>tol)
       c = (1+(eplas+e)/e0)^(1/n)*(e/(dt*edot0))^(1/m);
       f = sequiv/Y - 1.5*e*E/(Y*(1+nu))- c;
       dfde = -1.5*E/(Y*(1+nu)) - c*(1/(n*(eplas+e+e0)) + 1/(m*e));
       enew = e - f/dfde;
       if (enew<0)
%        e must be >0, so if new approx to e <0 the solution
%        must lie between current approx to e and zero.
         e = e/10;
       else
         e = enew;
       end
       err = abs(f);
     end
   end

end
