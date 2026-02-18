function [parentNode, parentQuadrant] = QTree_searchParentNode(childNode,QUADRANTS)
% PURPOSE: 
%  Search for the IDs of the parent cell / node per given quadrant
% INPUT: 
%  childNode = ID of child node
%  targetQuadrant  = which quadrant to search
% OUTPUT: 
%  parentNode
%  parentQuadrant
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

parentQuadrant = QUADRANTS(childNode);

switch parentQuadrant
    case 1
        parentNode = (childNode+2)/4;
    case 2
        parentNode = (childNode+1)/4;
    case 3
        parentNode = childNode/4;
    case 4
        parentNode = (childNode-1)/4;
end

parentQuadrant = QUADRANTS(parentNode);

return