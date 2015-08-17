function print_results(outfile,ndof,nnode,coords,nelem,connect,stress,eplas,dofs,nne)
%      Print nodal displacements, element strains && stresses:a file

fprintf(outfile,'Nodal Displacements: \n');
     fprintf(outfile,' Node      Coords         u1       u2 \n');
     for i = 1:nnode
      fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n', ...
                               i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
     end

   fprintf(outfile,'\n\n Strains and Stresses \n');

   lmncoord = zeros(ndof,nne);
   displacements = zeros(ndof,nne);

%
%   Loop over all the elements
%
   for lmn = 1:nelem

    fprintf(outfile,' \n Element; %d ',lmn);
    fprintf(outfile,'  \n int pt    Coords      eplas         e_11      e_22     e_12      s_11       s_22      s_12 \n');
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1:nne
        for i = 1:ndof
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
          displacements(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
      n = nne;
 
      npoints = 4;
      dNdx = zeros(n,ndof);
      dxdxi = zeros(ndof,ndof);
      strain = zeros(ndof,ndof);
      xi = zeros(ndof,1);
      x = zeros(ndof,1);
%
%  Set up integration points 
%
      xilist = integrationpoints(ndof,nne,npoints);
%
%  Loop over the integration points
%
     for intpt = 1:npoints

%     Compute shape functions && derivatives wrt local coords
%
       for i = 1:ndof
         xi(i) = xilist(i,intpt);
       end
       N = shapefunctions(n,ndof,xi);      
       dNdxi = shapefunctionderivs(n,ndof,xi);
%
%     Compute the coords of the integration point
%
      for i = 1:ndof
        x(i) = 0.;
        for a = 1:n
          x(i) = x(i) + lmncoord(i,a)*N(a);
        end
      end
%
%     Compute the jacobian matrix && its determinant
%
      for i = 1:ndof
        for j = 1:ndof
          dxdxi(i,j) = 0.;
          for a = 1:n
            dxdxi(i,j) = dxdxi(i,j) + lmncoord(i,a)*dNdxi(a,j);
          end
        end
      end

      dxidx = inv(dxdxi);
%
%     Convert shape function derivatives:derivatives wrt global coords
%
      for a = 1:n
        for i = 1:ndof
          dNdx(a,i) = 0.;
          for j = 1:ndof
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end
%
%     Compute the (infinitesimal) strain by differentiating displacements
%
      for i = 1:ndof
         for j = 1:ndof
            strain(i,j) = 0.;
            for a = 1:n
              strain(i,j) = strain(i,j) + 0.5*(displacements(i,a)*dNdx(a,j)+displacements(j,a)*dNdx(a,i));
            end
         end
      end

      fprintf(outfile,'%5d %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n', ...
        intpt,x(1),x(2),eplas(intpt,lmn),strain(1,1),strain(2,2),strain(1,2),stress(1,1,intpt,lmn),stress(2,2,intpt,lmn),stress(1,2,intpt,lmn));

     end
   end
end