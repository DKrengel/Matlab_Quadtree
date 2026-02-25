function [boxes1, boxes2, boxesLength] = QTree_searchNeighborBoxes(cell, ...
          nodeState, QUADRANTS, maxPairsCount, cellCenterX, cellCenterY, ...
          cellBottomX, cellBottomY, cellTopX, cellTopY, xBoxMax, yBoxMax)

% PURPOSE: 
%   Search neighboring cells of a given cell and collect neighboring bounding-box pairs.
% INPUT
%   cell             = ID of the reference cell
%   nodeState        = node state array (empty/internal/leaf-with-box-id)
%   QUADRANTS        = quadrant definition/lookup table used for traversal
%   maxPairsCount    = maximum number of bounding-box pairs to store
%   cellCenterX / cellCenterY
%                    = x / y coordinates of cell centers
%   cellBottomX / cellBottomY / cellTopX / cellTopY
%                    = spatial bounds (min x / min y / max x / max y) of each cell
%   xBoxMax / yBoxMax= maximum bounding-box width / height
% OUTPUT
%   boxes1 / boxes2  = indices of neighboring bounding-box pairs
%   boxesLength      = number of neighboring pairs found
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

isUncheckedCell = true(length(nodeState), 1);
isUncheckedCell(cell) = false;

MAX_RECHECK_CELLS = 100;
recheckCells = zeros(MAX_RECHECK_CELLS, 1);
minIndex = 0;
boxes1 = zeros(maxPairsCount, 1);
boxes2 = zeros(maxPairsCount, 1);
boxesLength = 0;

% First check: find all cells bordering the currently active one
adjacentCells = QTree_findAllAdjacentCells(cell, nodeState, QUADRANTS);
isUncheckedCell(adjacentCells) = false;

lengthAdjacentCells = length(adjacentCells);
BOX = nodeState(cell);

for i = 1:lengthAdjacentCells
    iCell = adjacentCells(i);
    boxTarget = nodeState(iCell);
    if boxTarget > 0 % if b is leaf node
        boxesLength = boxesLength + 1;
        if boxTarget < BOX
            boxes1(boxesLength) = boxTarget;
            boxes2(boxesLength) = BOX;
        else
            boxes1(boxesLength) = BOX;
            boxes2(boxesLength) = boxTarget;
        end
    end
end
maxIndex = lengthAdjacentCells;
recheckCells(1:lengthAdjacentCells) = adjacentCells;

while maxIndex > minIndex
    minIndex = minIndex + 1;
    iCell = recheckCells(minIndex);
    cells = QTree_findAllAdjacentCells(iCell, nodeState, QUADRANTS);
    for j = 1:length(cells)
        jCell = cells(j);
        boxTarget = nodeState(jCell);
        if isUncheckedCell(jCell)
            if boxTarget > 0 % if b is leaf node
                boxesLength = boxesLength + 1;
                if boxTarget < BOX
                    boxes1(boxesLength) = boxTarget;
                    boxes2(boxesLength) = BOX;
                else
                    boxes1(boxesLength) = BOX;
                    boxes2(boxesLength) = boxTarget;
                end
            else

                baseCell = cell;
                targetCell = jCell;

                if cellCenterX(baseCell) < cellCenterX(targetCell)
                    xBase = cellTopX(baseCell);
                    xTarget = cellBottomX(targetCell);
                    isNeighborX = xBoxMax > (xTarget - xBase);
                else
                    xBase = cellBottomX(baseCell);
                    xTarget = cellTopX(targetCell);
                    isNeighborX = xBoxMax > (xBase - xTarget);
                end

                if cellCenterY(baseCell) < cellCenterY(targetCell)
                    yBase = cellTopY(baseCell);
                    yTarget = cellBottomY(targetCell);
                    isNeighborY = yBoxMax > (yTarget - yBase);
                else
                    yBase = cellBottomY(baseCell);
                    yTarget = cellTopY(targetCell);
                    isNeighborY = yBoxMax > (yBase - yTarget);
                end

                isNeighborCell = false;
                if isNeighborX && isNeighborY
                    isNeighborCell = true;
                end


                if isNeighborCell
                    % Close cell-to-cell distance
                    maxIndex = maxIndex + 1;
                    recheckCells(maxIndex) = jCell;  
                end
            end
            isUncheckedCell(jCell) = false;
        end
    end
end

boxes1 = boxes1(1:boxesLength);
boxes2 = boxes2(1:boxesLength);

return