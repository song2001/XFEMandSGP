function [XY,XY1,LE,LE1] = ABAQUSmesh2D(Q)
%
% Read from ABAQUS generated mesh (2D version)
%
% Emiio Martínez Pañeda (martinezemilio@uniovi.es)
%
%
%--------------------------------------------------------------------------

% Nodal coordiantes (nodes.txt is obtained from ABAQUS input file by
% extracting the text after *Node and suppressing the commas)
fileID = fopen('nodes.txt','r');
formatSpec = '%d %f %f';
sizeA = [3 Inf];
Ab = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
XY=Ab';
XY(:,1)=[];

% Element connectivity matrix (conec.txt is obtained from ABAQUS input file
% by extracting the next after *Element and suppressing the commas)

fileID = fopen('conec.txt','r');
if Q(1)==4 && Q(2)==4 % Linear FEM and Linear X-FEM
formatSpec = '%d %d %d %d %d %d %d %d %d';
sizeA = [9 Inf];
Ac = fscanf(fileID,formatSpec,sizeA);
LE=Ac';
LE(:,1)=[];
LE(:,5)=[];LE(:,5)=[];LE(:,5)=[];LE(:,5)=[];
XY=XY(1:max(LE(:)),:);
LE1=LE;XY1=XY;
elseif Q(1)==8 && Q(2)==4 % Quadratic FEM and Linear X-FEM
formatSpec = '%d %d %d %d %d %d %d %d %d';
sizeA = [9 Inf];
Ac = fscanf(fileID,formatSpec,sizeA);
LE=Ac';
LE(:,1)=[];
LE1=LE;
LE1(:,5)=[];LE1(:,5)=[];LE1(:,5)=[];LE1(:,5)=[];
XY1=XY(1:max(LE1(:)),:);
elseif Q(1)==8 && Q(2)==8 % Quadratic FEM and Quadratic X-FEM    
formatSpec = '%d %d %d %d %d %d %d %d %d';
sizeA = [9 Inf];
Ac = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
LE=Ac';
LE(:,1)=[];
LE1=LE;XY1=XY;
end

end