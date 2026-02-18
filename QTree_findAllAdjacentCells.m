function adjacentCells = QTree_findAllAdjacentCells(baseCell,nodeState,QUADRANTS)
% PURPOSE: 
%  find the IDs of the 8 cells bordering to the current one
% INPUT: 
%  baseCell  = which cell is currently active as reference point for neighbors
%  nodeState
%  QUADRANTS
% OUTPUT: 
%  adjacentCells
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

MAX_ADJACENT_CELLS = 100;
adjacentCells = zeros(1, MAX_ADJACENT_CELLS);
count = 1;

for direction = 1:4
    %4 edge contacts (up/down/left/right)
    edgeNeighborCells = QTree_findEdgeAdjacentCell(baseCell,direction,nodeState,QUADRANTS);
    if edgeNeighborCells(1) ~= 0 % NOTE: If there is an edge-cell
        adjacentCells(count:(count+length(edgeNeighborCells)-1)) = edgeNeighborCells;
        count = count + length(edgeNeighborCells);
    end
    %4 corner contacts (left up/ right up/ right down/ left down)
    cornerAdjacentCell = QTree_findCornerAdjacentCell(baseCell,direction,nodeState,QUADRANTS);

    if cornerAdjacentCell ~= 0
        adjacentCells(count) = cornerAdjacentCell;
        count = count + 1;
    end
end

adjacentCells = nonzeros(adjacentCells).';

return