function [objpos,volsize]=findobject(image,connection)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Sébastien Huet
%   Email: sebastien.huet (at) univ-rennes1.fr
%   Date: 29-07-2021
%   
%   This function is used in the main routines Recruitment_auto3 and PATag_decondense_3 to collect the position and sizes of the segmented objects
%
%%%%%%%%%%%%%%%%%%%%%%%%%


[objects,objnum]=bwlabeln(image,connection);
objpos=cell(objnum,1);
volsize=zeros(objnum,1);
for i=1:objnum
    objpos{i}=find(objects==i);
    volsize(i)=size(objpos{i},1);
end