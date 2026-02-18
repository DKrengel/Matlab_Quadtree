function edgeNeighborCells = QTree_findDeeperEdgeCells(baseEdgeAdjacentCell, edgeDirection, nodeState)

% PURPOSE: 
% INPUT: 
%  baseEdgeAdjacentCell 
%  edgeDirection
%  nodeState
% OUTPUT: 
%  edgeNeighborCells
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