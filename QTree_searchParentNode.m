function [parentNode, parentQuadrant] = QTree_searchParentNode(childNode,QUADRANTS)
% PURPOSE: 
%   Find the parent node of a given child node and determine its quadrant position.
% INPUT
%   childNode   = ID of the child node
%   QUADRANTS   = quadrant definition/lookup table used for traversal
% OUTPUT
%   parentNode      = ID of the parent node
%   parentQuadrant  = quadrant index of the child node within its parent
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