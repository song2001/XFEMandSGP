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

botNodes = find(XY(:,2)==miny);
rightNodes = find(XY(:,1)==maxx);
topNodes = find(XY(:,2)==maxy);
leftNodes = find(XY(:,1)==minx);


if NNE==4
p=1; store=zeros(NE,1);    
%rightEdge(m1,:) = [i,LE(j,n)] ;
 for i=1:size(botNodes)
 [r,~]=find(LE == botNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store(r(1))==0
   store(r(1))=botNodes(i);
   else
   botEdge(p,:)=[store(r(1)) botNodes(i)];
   p=p+1;
   end
   if store(r(2))==0
   store(r(2))=botNodes(i);
   else
   botEdge(p,:)=[store(r(2)) botNodes(i)];
   p=p+1;
   end
  else %corner node
   if store(r)==0
   store(r)=botNodes(i);
   else
   botEdge(p,:)=[store(r) botNodes(i)];
   p=p+1;
   end   
  end
 end
 
p=1; store=zeros(NE,1);
 for i=1:size(rightNodes)
 [r,~]=find(LE == rightNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store(r(1))==0
   store(r(1))=rightNodes(i);
   else
   rightEdge(p,:)=[store(r(1)) rightNodes(i)];
   p=p+1;
   end
   if store(r(2))==0
   store(r(2))=rightNodes(i);
   else
   rightEdge(p,:)=[store(r(2)) rightNodes(i)];
   p=p+1;
   end
  else %corner node
   if store(r)==0
   store(r)=rightNodes(i);
   else
   rightEdge(p,:)=[store(r) rightNodes(i)];
   p=p+1;
   end   
  end
 end

p=1; store=zeros(NE,1);
 for i=1:size(topNodes)
 [r,~]=find(LE == topNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store(r(1))==0
   store(r(1))=topNodes(i);
   else
   topEdge(p,:)=[store(r(1)) topNodes(i)];
   p=p+1;
   end
   if store(r(2))==0
   store(r(2))=topNodes(i);
   else
   topEdge(p,:)=[store(r(2)) topNodes(i)];
   p=p+1;
   end
  else %corner node
   if store(r)==0
   store(r)=topNodes(i);
   else
   topEdge(p,:)=[store(r) topNodes(i)];
   p=p+1;
   end   
  end
 end 

p=1; store=zeros(NE,1);
 for i=1:size(leftNodes)
 [r,~]=find(LE == leftNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store(r(1))==0
   store(r(1))=leftNodes(i);
   else
   leftEdge(p,:)=[store(r(1)) leftNodes(i)];
   p=p+1;
   end
   if store(r(2))==0
   store(r(2))=leftNodes(i);
   else
   leftEdge(p,:)=[store(r(2)) leftNodes(i)];
   p=p+1;
   end
  else %corner node
   if store(r)==0
   store(r)=topNodes(i);
   else
   leftEdge(p,:)=[store(r) leftNodes(i)];
   p=p+1;
   end   
  end
 end  
elseif NNE==8
p=1; store1=zeros(NE,1); store2=zeros(NE,1);    
 for i=1:size(botNodes)
 [r,~]=find(LE == botNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store1(r(1))==0
   store1(r(1))=botNodes(i);
   elseif store2(r(1))==0
   store2(r(1))=botNodes(i);   
   else
   botEdge(p,:)=[store1(r(1)) store2(r(1)) botNodes(i)];
   p=p+1;
   end
   if store1(r(2))==0
   store1(r(2))=botNodes(i);
   elseif store2(r(2))==0
   store2(r(2))=botNodes(i);   
   else
   botEdge(p,:)=[store1(r(2)) store2(r(2)) botNodes(i)];
   p=p+1;
   end
  else
   if store1(r)==0
   store1(r)=botNodes(i);
   elseif store2(r)==0
   store2(r)=botNodes(i);    
   else    
   botEdge(p,:)=[store1(r) store2(r) botNodes(i)];
   p=p+1;
   end   
  end
 end
 p=1; store1=zeros(NE,1); store2=zeros(NE,1);
 for i=1:size(rightNodes)
 [r,~]=find(LE == rightNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store1(r(1))==0
   store1(r(1))=rightNodes(i);
   elseif store2(r(1))==0
   store2(r(1))=rightNodes(i);   
   else
   rightEdge(p,:)=[store1(r(1)) store2(r(1)) rightNodes(i)];
   p=p+1;
   end
   if store1(r(2))==0
   store1(r(2))=rightNodes(i);
   elseif store2(r(2))==0
   store2(r(2))=rightNodes(i);   
   else
   rightEdge(p,:)=[store1(r(2)) store2(r(2)) rightNodes(i)];
   p=p+1;
   end
  else
   if store1(r)==0
   store1(r)=rightNodes(i);
   elseif store2(r)==0
   store2(r)=rightNodes(i);    
   else    
   rightEdge(p,:)=[store1(r) store2(r) rightNodes(i)];
   p=p+1;
   end   
  end
 end

 p=1; store1=zeros(NE,1); store2=zeros(NE,1);
 for i=1:size(topNodes)
 [r,~]=find(LE == topNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store1(r(1))==0
   store1(r(1))=topNodes(i);
   elseif store2(r(1))==0
   store2(r(1))=topNodes(i);   
   else
   topEdge(p,:)=[store1(r(1)) store2(r(1)) topNodes(i)];
   p=p+1;
   end
   if store1(r(2))==0
   store1(r(2))=topNodes(i);
   elseif store2(r(2))==0
   store2(r(2))=topNodes(i);   
   else
   topEdge(p,:)=[store1(r(2)) store2(r(2)) topNodes(i)];
   p=p+1;
   end
  else
   if store1(r)==0
   store1(r)=topNodes(i);
   elseif store2(r)==0
   store2(r)=topNodes(i);    
   else    
   topEdge(p,:)=[store1(r) store2(r) topNodes(i)];
   p=p+1;
   end   
  end
 end 
 
 p=1; store1=zeros(NE,1); store2=zeros(NE,1);
 for i=1:size(leftNodes)
 [r,~]=find(LE == leftNodes(i));
 [c, ~] = size(r);
  if c==2 
   if store1(r(1))==0
   store1(r(1))=leftNodes(i);
   elseif store2(r(1))==0
   store2(r(1))=leftNodes(i);   
   else
   leftEdge(p,:)=[store1(r(1)) store2(r(1)) leftNodes(i)];
   p=p+1;
   end
   if store1(r(2))==0
   store1(r(2))=leftNodes(i);
   elseif store2(r(2))==0
   store2(r(2))=leftNodes(i);   
   else
   leftEdge(p,:)=[store1(r(2)) store2(r(2)) leftNodes(i)];
   p=p+1;
   end
  else
   if store1(r)==0
   store1(r)=leftNodes(i);
   elseif store2(r)==0
   store2(r)=leftNodes(i);    
   else    
   leftEdge(p,:)=[store1(r) store2(r) leftNodes(i)];
   p=p+1;
   end   
  end
 end  
end

bnodes = {botNodes rightNodes topNodes leftNodes} ;
enodes = {botEdge rightEdge topEdge leftEdge} ;
end