function edgeAdjacentCells = QTree_findEdgeAdjacentCell(edgeAdjacentCell,direction,nodeState,QUADRANTS)

% PURPOSE: 
%  find the IDs of the cell contacting the edges of the currently active cell in Quadtrant
% INPUT: 
%  edgeAdjacentCell (singular)
%  direction
%  nodeState
%  QUADRANTS
% OUTPUT: 
%  edgeAdjacentCells (plural!!)
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

REVERSE_EDGE_DIRECTION = [3;4;1;2];
IS_NOT_ADJACENT = [
    1 1 0 0;
    0 1 0 1;
    0 0 1 1;
    1 0 1 0
];
REFLECT = [
    0  0 -2 -2;
    1  0  1  0;
    2  2  0  0;
    0 -1  0 -1
];
FACE = [
    3 4 1 2;
    2 1 4 3;
    3 4 1 2;
    2 1 4 3
];

up = 0;
neighborQuadrant = QUADRANTS(edgeAdjacentCell); % node number
isSearchNeeded   = true;
targetQuadrants  = zeros(10,1);

while IS_NOT_ADJACENT(direction,neighborQuadrant)
    symmetricalQuadrant = neighborQuadrant;
    [edgeAdjacentCell,neighborQuadrant] = QTree_searchParentNode(edgeAdjacentCell,QUADRANTS);

    if edgeAdjacentCell == 1
        edgeAdjacentCells = 0; % NOTE: Neighborhood node does not exist.
        isSearchNeeded = false;
        break
    end
    up = up + 1;
    targetQuadrants(up) = FACE(direction,symmetricalQuadrant);
end

if isSearchNeeded
    isDownNeeded = true;
    edgeAdjacentCell = edgeAdjacentCell + REFLECT(direction,neighborQuadrant);
    if nodeState(edgeAdjacentCell) >= 0 % NOTE: Empty-node or Leaf-node
        edgeAdjacentCells = edgeAdjacentCell;
        isDownNeeded = false;
    end

    % while isDownNeeded & (is_cor_neighbor==false)
    while isDownNeeded
        nodeFeature = nodeState(edgeAdjacentCell);

        if up <= 0
            if nodeFeature == -1
                edgeAdjacentCells = QTree_findDeeperEdgeCells(edgeAdjacentCell, REVERSE_EDGE_DIRECTION(direction), nodeState);
            else
                edgeAdjacentCells = edgeAdjacentCell;
            end

            break
        end

        if nodeFeature >= 0 % NOTE: Empty-node or Leaf-node
            edgeAdjacentCells = edgeAdjacentCell;
            break
        else
            edgeAdjacentCell = QTree_searchChildNode(edgeAdjacentCell,targetQuadrants(up));
        end

        up = up - 1;
    end
end

return