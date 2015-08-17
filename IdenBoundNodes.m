function [bnodes,enodes] = IdenBoundNodes(XY,LE)
% Get the number of elements in the X- and Y- directions
% (Note that we will always work with structured meshes)
x=XY(LE(1,1),1);
y=XY(LE(1,1),2);
nnx=-1; nny=-1; 

maxx=XY(LE(1,1),1);
maxy=XY(LE(1,1),2);
minx=maxx; miny=maxy;
[NUM, ~] = size(XY);
[NE, NNE] = size(LE);
for i=1:NUM
 if XY(i,1)>maxx  
  maxx=XY(i,1);
 end
 if XY(i,1)<minx
  minx=XY(i,1);
 end 
 if XY(i,2)>maxy  
  maxy=XY(i,2);
 end 
 if XY(i,2)<miny
  miny=XY(i,2);
 end
end 
m1=1;m2=1;m3=1;m4=1;
if NNE==4
 for i=1:NUM
     
  if XY(i,1)==maxx
   for j=1:NE
    for k=1:4   
     if LE(j,k)==i
      for n=1:4
       if ((XY(LE(j,n),1)==maxx) && (n~=k))
       rightEdge(m1,:) = [i,LE(j,n)] ;
       m1=m1+1;
       end  
      end
     end
    end
   end    
  end
  
  if XY(i,1)==minx
   for j=1:NE
    for k=1:4   
     if LE(j,k)==i
      for n=1:4
       if ((XY(LE(j,n),1)==minx) && (n~=k))
       leftEdge(m2,:) = [i,LE(j,n)] ;
       m2=m2+1;
       end  
      end
     end
    end
   end    
  end  
  
  if XY(i,2)==maxy
   for j=1:NE
    for k=1:4   
     if LE(j,k)==i
      for n=1:4
       if ((XY(LE(j,n),2)==maxy) && (n~=k))
       topEdge(m3,:) = [i,LE(j,n)] ;
       m3=m3+1;
       end  
      end
     end
    end
   end    
  end
  
  if XY(i,2)==miny
   for j=1:NE
    for k=1:4   
     if LE(j,k)==i
      for n=1:4
       if ((XY(LE(j,n),2)==miny) && (n~=k))
       botEdge(m4,:) = [i,LE(j,n)] ;
       m4=m4+1;
       end  
      end
     end
    end
   end    
  end    
  
  
 end  
end

rightEdge(2:2:end,:) = [];
leftEdge(2:2:end,:) = [];
topEdge(2:2:end,:) = [];
botEdge(2:2:end,:) = [];

botNodes = unique(botEdge) ;
rightNodes = unique(rightEdge) ;
topNodes = unique(topEdge) ;
leftNodes = unique(leftEdge) ;

bnodes = {botNodes rightNodes topNodes leftNodes} ;
enodes = {botEdge rightEdge topEdge leftEdge} ;
end