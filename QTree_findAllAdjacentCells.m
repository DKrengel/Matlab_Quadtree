function adjacentCells = QTree_findAllAdjacentCells(baseCell,nodeState,QUADRANTS)
% PURPOSE: 
%  find the IDs of the all cells bordering to the current one in the 8 directions 
%  (N, S, E, W, NE, NW, SE, SW) 
% INPUT: 
%  baseCell  = which cell is currently active as reference point for neighbors
%  nodeState = node state array (empty/internal/leaf-with-box-id)
%  QUADRANTS = quadrant definition/lookup table used for traversal
% OUTPUT: 
%  adjacentCells = list of IDs of adjacent cells (8 directions)
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