function isNeighborCell = QTree_isNeighborCell(baseCell, targetCell, ...
           cellCenterX, cellCenterY, cellBottomX, cellBottomY, ...
           cellTopX, cellTopY, xBoxMax, yBoxMax)

% PURPOSE
%   Check whether a target cell is spatially adjacent to a base cell.
% INPUT
%   baseCell    = ID of the reference cell
%   targetCell  = ID of the candidate neighboring cell
%   cellCenterX / cellCenterY
%               = x / y coordinates of cell centers
%   cellBottomX / cellBottomY / cellTopX / cellTopY
%               = spatial bounds (min x / min y / max x / max y) of each cell
%   xBoxMax / yBoxMax
%               = maximum bounding-box width / height
% OUTPUT
%   isNeighborCell = logical value indicating whether the target cell is adjacent to the base cell
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

if cellCenterX(baseCell) < cellCenterX(targetCell)
    xBase   = cellTopX(baseCell);
    xTarget = cellBottomX(targetCell);
    isNeighborX = xBoxMax > (xTarget - xBase);
else
    xBase   = cellBottomX(baseCell);
    xTarget = cellTopX(targetCell);
    isNeighborX = xBoxMax > (xBase - xTarget);
end

if cellCenterY(baseCell) < cellCenterY(targetCell)
    yBase   = cellTopY(baseCell);
    yTarget = cellBottomY(targetCell);
    isNeighborY = yBoxMax > (yTarget - yBase);
else
    yBase   = cellBottomY(baseCell);
    yTarget = cellTopY(targetCell);
    isNeighborY = yBoxMax > (yBase - yTarget);
end

isNeighborCell = false;
if isNeighborX && isNeighborY
    isNeighborCell = true;
end

return