function edgeNeighborCells = QTree_findDeeperEdgeCells(baseEdgeAdjacentCell, edgeDirection, nodeState)
% PURPOSE:
%   Find deeper edge-neighbor cells of a given edge-adjacent cell in the specified direction.
% INPUT: 
%   baseEdgeAdjacentCell = ID of the base edge-adjacent cell
%   edgeDirection        = edge direction identifier (e.g., N, S, E, W)
%   nodeState            = node state array (empty/internal/leaf-with-box-id)
% OUTPUT: 
%   edgeNeighborCells    = list of IDs of deeper edge-neighbor cells
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

MAX_INTERNAL_NODES = 100;
MAX_NEIGHBOR_NODES = 100;

internalNodes     = zeros(MAX_INTERNAL_NODES, 1);
edgeNeighborCells = zeros(MAX_NEIGHBOR_NODES, 1);

internalNodes(1)   = baseEdgeAdjacentCell;
internalNodesIndex = 1;
internalNodesIndexEnd = 1;
neighborNodesIndex    = 1;

while internalNodes(internalNodesIndex) > 0
    parentCell = internalNodes(internalNodesIndex);
    childCells = QTree_childCellsEachDirection(parentCell, edgeDirection);
    for childCell = childCells
        if nodeState(childCell) == -1
            internalNodesIndexEnd = internalNodesIndexEnd + 1;
            internalNodes(internalNodesIndexEnd) = childCell;
        else
            edgeNeighborCells(neighborNodesIndex) = childCell;
            neighborNodesIndex = neighborNodesIndex + 1;
        end
    end
    internalNodesIndex = internalNodesIndex + 1;
end

neighborNodesCount = neighborNodesIndex - 1;

edgeNeighborCells = edgeNeighborCells(1:neighborNodesCount);

return