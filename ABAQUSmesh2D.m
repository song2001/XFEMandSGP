function [XY,LE] = ABAQUSmesh2D(Q)
%
% Read from ABAQUS generated mesh (2D version - CPE8R element)
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
if Q==8
formatSpec = '%d %d %d %d %d %d %d %d %d';
sizeA = [9 Inf];
else
formatSpec = '%d %d %d %d %d';
sizeA = [5 Inf];
end
Ac = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
LE=Ac';
LE(:,1)=[];

end