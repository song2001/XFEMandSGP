function list = facenodes(nelnodes,face)

  i4 = [2,3,4,1]; 

  list = zeros(nfacenodes(nelnodes),1);
  
     if (nelnodes==4) 
       list(1) = face;
       list(2) = i4(face);
     elseif (nelnodes==8) 
       list(1) = face;
       list(2) = i4(face);
       list(3) = face+4;
     end
   
end