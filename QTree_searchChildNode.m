function childNode = QTree_searchChildNode(parentNode,targetQuadrant)
% PURPOSE: 
%  Search for the IDs of the child cells / nodes per given quadrant
% INPUT: 
%  parentNode = ID of parent node
%  targetQuadrant  = which quadrant to search
% OUTPUT: 
%  childNode
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

switch targetQuadrant
    case 1
        childNode = 4*parentNode-2;
    case 2
        childNode = 4*parentNode-1;
    case 3
        childNode = 4*parentNode;
    case 4
        childNode = 4*parentNode+1;
    otherwise
        childNode = 0;
end

return
