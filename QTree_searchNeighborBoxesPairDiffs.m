function [boxes1, boxes2, boxesLength] = QTree_searchNeighborBoxesPairDiffs(cell, nodeState, QUADRANTS, maxPairsCount,...
    cellCenterX, cellCenterY, cellBottomX, cellBottomY, cellTopX, cellTopY, xBoxMax, yBoxMax, fi_n, la_n)

% PURPOSE: 
% Neighborhood computation via a quadtree-approach
% INPUT: 
% OUTPUT: 
%  boxes1, boxes2 = 
%  boxesLength = 
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

GRADIENT_22_5 = sqrt(2) - 1;
GRADIENT_67_5 = sqrt(2) + 1;

FIRST_NODE_FLIP = flip(fi_n);
LAST_NODE_FLIP = flip(la_n);

LENGTH_NODE_STATE = length(nodeState);

% MEMO: Scope expands, but defining it here makes it faster
isUncheckedCell = true(LENGTH_NODE_STATE, 1);
isUncheckedCell(cell) = false;

% First check
[adjacentCells] = QTree_findAllAdjacentCells(cell, nodeState, QUADRANTS);

MAX_RECHECK_CELLS = 100;
recheckCells = zeros(MAX_RECHECK_CELLS, 1);
minIndex = 0;
boxes1 = zeros(maxPairsCount, 1);
boxes2 = zeros(maxPairsCount, 1);
boxesLength = 0;

isUncheckedCell(adjacentCells) = false;

lengthAdjacentCells = length(adjacentCells);
BOX = nodeState(cell);

% If the 9 cells centered on the base-cell(cell variable)
% have the same size (same depth), the search area can be narrowed
depthFlip       = find((FIRST_NODE_FLIP <= cell), 1);
areSameDepth    = FIRST_NODE_FLIP(depthFlip) <= adjacentCells & adjacentCells <= LAST_NODE_FLIP(depthFlip);
areAllSameDepth = all(areSameDepth);

for i = 1:lengthAdjacentCells
    iCell     = adjacentCells(i);
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
    iCell    = recheckCells(minIndex);

    % START: calcNodeDirection
    % input: node
    % output: directions

    baseX = cellCenterX(cell);
    baseY = cellCenterY(cell);
    targetX = cellCenterX(iCell);
    targetY = cellCenterY(iCell);

    gradient = (targetY - baseY) / (targetX - baseX);

    % Orientation: right and left
    if gradient <= GRADIENT_22_5 && -GRADIENT_22_5 <= gradient
        if baseX < targetX
            if areAllSameDepth
                directions = [2 2 4];
                isCorner = [true false true];
            else
                directions = [1 2 2 4 3];
                isCorner = [false true false true false];
            end
        else
            if areAllSameDepth
                directions = [3 4 1];
                isCorner = [true false true];
            else
                directions = [3 3 4 1 1];
                isCorner = [false true false true false];
            end
        end
    % Orientation: top-right and bottom-left
    elseif GRADIENT_22_5 < gradient && gradient < GRADIENT_67_5
        if baseX < targetX
            if areAllSameDepth
                directions = [1 2 2];
                isCorner = [false true false];
            else
                directions = [1 1 2 2 4];
                isCorner = [true false true false true];
            end
        else
            if areAllSameDepth
                directions = [3 3 4];
                isCorner = [false true false];
            else
                directions = [4 3 3 4 1];
                isCorner = [true false true false true];
            end
        end
    % Orientation: bottom-right and top-left
    elseif -GRADIENT_67_5 < gradient && gradient < -GRADIENT_22_5
        if baseX < targetX
            if areAllSameDepth
                directions = [2 4 3];
                isCorner = [false true false];
            else
                directions = [2 2 4 3 3];
                isCorner = [true false true false true];
            end
        else
            if areAllSameDepth
                directions = [4 1 1];
                isCorner = [false true false];
            else
                directions = [3 4 1 1 2];
                isCorner = [true false true false true];
            end
        end
    % Orientation: top and bottom
    else
        if baseY < targetY
            if areAllSameDepth
                directions = [1 1 2];
                isCorner = [true false true];
            else
                directions = [4 1 1 2 2];
                isCorner = [false true false true false];
            end
        else
            if areAllSameDepth
                directions = [4 3 3];
                isCorner = [true false true];
            else
                directions = [2 4 3 3 4];
                isCorner = [false true false true false];
            end
        end
    end

    countDirections = length(directions);

    % End: calcNodeDirection

    % START: calcNeighborNodesByDirection
    % input: node, directions
    % output: nodes

    directedNeighborNodes = zeros(20, 1);
    directedNeighborNodesIndex = 0;

    for d = 1:countDirections
        if isCorner(d)
            direction = directions(d);

            cornerAdjacentCell = QTree_findCornerAdjacentCell(iCell, direction, nodeState, QUADRANTS);

            if cornerAdjacentCell ~= 0
                directedNeighborNodesIndex = directedNeighborNodesIndex + 1;
                directedNeighborNodes(directedNeighborNodesIndex) = cornerAdjacentCell;
            end
        else
            direction = directions(d);
            edgeNeighborCells = QTree_findEdgeAdjacentCell(iCell, direction, nodeState, QUADRANTS);

            for j = 1:length(edgeNeighborCells)
                edgeNeighborCell = edgeNeighborCells(j);
                if edgeNeighborCell ~= 0
                    directedNeighborNodesIndex = directedNeighborNodesIndex + 1;
                    directedNeighborNodes(directedNeighborNodesIndex) = edgeNeighborCell;
                end
            end
        end
    end

    % END: calcNeighborNodesByDirection

    for j = 1:directedNeighborNodesIndex
        jCell = directedNeighborNodes(j);
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
                isNeighbor = QTree_isNeighborCell(cell, jCell, cellCenterX, cellCenterY,...
                          cellBottomX, cellBottomY, cellTopX, cellTopY, xBoxMax, yBoxMax);
                if isNeighbor
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