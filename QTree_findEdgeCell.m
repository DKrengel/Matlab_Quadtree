function [cornerAdjacentCell,hasPolygonInEdge] = QTree_findEdgeCell(cornerAdjacentCell,edgeDirection,nodeState,QUADRANTS)

% PURPOSE: 
% INPUT: 
%  cornerAdjacentCell
%  edgeDirection
%  nodeState
%  QUADRANTS
% OUTPUT: 
%  cornerAdjacentCell
%  hasPolygonInEdge
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

% For calculating quadrant tangent to the specified direction
REFLECT = [
    0  0 -2 -2;
    1  0  1  0;
    2  2  0  0;
    0 -1  0 -1
];
% Quadrant tangent to the specified direction
FACE = [
    3 4 1 2;
    2 1 4 3;
    3 4 1 2;
    2 1 4 3
];
% TODO: Write comments!!!
IS_NOT_ADJACENT = [
    1 1 0 0;
    0 1 0 1;
    0 0 1 1;
    1 0 1 0
];

up = 0;
neighborQuadrant = QUADRANTS(cornerAdjacentCell); % node number
isSearchNeeded   = true;
targetQuadrants  = zeros(10,1);
hasPolygonInEdge = false;

while IS_NOT_ADJACENT(edgeDirection,neighborQuadrant)
    symmetricalQuadrant = neighborQuadrant;
    [cornerAdjacentCell,neighborQuadrant] = QTree_searchParentNode(cornerAdjacentCell,QUADRANTS);

    if cornerAdjacentCell == 1
        cornerAdjacentCell = 0; % NOTE: Neighborhood node does not exist.
        isSearchNeeded = false;
        break
    end
    up = up + 1;
    targetQuadrants(up) = FACE(edgeDirection,symmetricalQuadrant);
end

if isSearchNeeded
    isDownNeeded = true;
    cornerAdjacentCell = cornerAdjacentCell + REFLECT(edgeDirection,neighborQuadrant);
    switch nodeState(cornerAdjacentCell)
        case 0 % if empty
            isDownNeeded = false;
        case -1
            % Not implemented
        otherwise
            isDownNeeded = false;
            hasPolygonInEdge = true;
    end

    while isDownNeeded
        nodeFeature = nodeState(cornerAdjacentCell);
        if up <= 0
            if nodeFeature ~= 0 && nodeFeature ~= -1
                hasPolygonInEdge = true;
            end
            break
        end
        switch nodeFeature
            case 0
                break
            case -1
                cornerAdjacentCell = QTree_searchChildNode(cornerAdjacentCell,targetQuadrants(up));
            otherwise
                hasPolygonInEdge = true;
                break
        end
        up = up - 1;
    end
end

return