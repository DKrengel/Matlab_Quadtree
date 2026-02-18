function isNeighborCell = QTree_isNeighborCell(baseCell, targetCell, ...
           cellCenterX, cellCenterY, cellBottomX, cellBottomY, ...
           cellTopX, cellTopY, xBoxMax, yBoxMax)

% PURPOSE: 
% INPUT: 
%  baseCell, 
%  targetCell
%  cellCenterX, cellCenterY
%  cellBottomX, cellBottomY
%  cellTopX, cellTopY
%  xBoxMax, yBoxMax
% OUTPUT: 
%  isNeighborCell
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