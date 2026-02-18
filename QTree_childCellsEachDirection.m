function childCells = childCellsEachDirection(parentCell, edgeDirection)

% PURPOSE: 
% INPUT: 
%  parentCell 
%  edgeDirection
% OUTPUT: 
%  childCells
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori


switch edgeDirection
    case 1
        childCells = [4*parentCell-2 4*parentCell-1];
    case 2
        childCells = [4*parentCell-1 4*parentCell+1];
    case 3
        childCells = [4*parentCell 4*parentCell+1];
    case 4
        childCells = [4*parentCell-2 4*parentCell];
end

return