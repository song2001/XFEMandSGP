% Function to perform Voigt mapping algorithm
function [Tensor]=Voigtmap(Voigt)
dl=[[1,0,0];[0,1,0];[0,0,1]];
Tensor=zeros(3,3,3,3);
     for i = 1 : 3
       for j = 1 : 3
         for k = 1 : 3
           for l = 1 : 3
             p=i*dl(i,j)+(1-dl(i,j))*(9-i-j);
             q=k*dl(k,l)+(1-dl(k,l))*(9-k-l);
             Tensor(i,j,k,l) = Voigt(p,q);
            end
          end
        end
     end
end