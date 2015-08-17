PROP=[260000, 0.3];

% Elasticity (4th order)

C = zeros(2,2,2,2);
  
     for i = 1:2
       for j = 1:2
         for k = 1:2
           for l = 1:2 
             if (i==j && k==l)
                   C(i,j,k,l) = C(i,j,k,l)+PROP(2)*PROP(1)/((1+PROP(2))*(1-2*PROP(2)));
             end
             if (i==l && k==j)
                 C(i,j,k,l) = C(i,j,k,l)+PROP(1)/(2*(1+PROP(2)));
             end
             if (i==k && j==l)
                 C(i,j,k,l) = C(i,j,k,l)+PROP(1)/(2*(1+PROP(2)));
             end
           end
         end
       end
     end
     
% Elasticity (2nd order)     

  LAM=PROP(1)*PROP(2)/((1+PROP(2))*(1-2*PROP(2)));
  MU=PROP(1)/(2*(1+PROP(2)));
  %
  %
  ETAN=[LAM+2*MU LAM       0;
        LAM      LAM+2*MU  0;
        0        0        MU];   
    
% Elasticity 4th order from 2nd order

C1= zeros(2,2,2,2);
C1(1,1,1,1)=ETAN(1,1);
C1(2,2,2,2)=ETAN(2,2);
C1(1,1,2,2)=ETAN(1,2);
C1(2,2,1,1)=ETAN(2,1);
C1(1,2,1,2)=ETAN(3,3);
C1(2,1,1,2)=ETAN(3,3);
C1(1,2,2,1)=ETAN(3,3);
C1(2,1,2,1)=ETAN(3,3);
     for i = 1:2
       for j = 1:2
         for k = 1:2
           for l = 1:2 
            if (C(i,j,k,l) ~=  C1(i,j,k,l))
             fprintf(1,'\n Error \n');
            end    
           end
         end
       end
     end
