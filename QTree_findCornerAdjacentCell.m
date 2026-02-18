function cornerAdjacentCell = QTree_findCornerAdjacentCell(cornerAdjacentCell,direction,nodeState,QUADRANTS)

% PURPOSE: 
%  find the IDs of the cell contacting the corners of the currently active cell in Quadtrant
% INPUT: 
%  cornerAdjacentCell
%  direction
%  nodeState
%  QUADRANTS
% OUTPUT: 
%  cornerAdjacentCell
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

HAS_SAME_PARENT = [
    % How to use: HAS_SAME_PARENT(quadrant, direction)
    0 0 0 1;
    0 0 1 0;
    0 1 0 0;
    1 0 0 0
];

DELTA_QUADRANT_HAS_SAME_PARENT = [3 1 -1 -3];
PARENT_NODE_ADJACENT_DIRECTION = [
    % How to use: PARENT_NODE_ADJACENT_DIRECTION(quadrant, direction)
    0 1 4 0;
    1 0 0 2;
    4 0 0 3;
    0 2 3 0
];

SYMMETRIC_CORNER_DIRECTION = [4;3;2;1];

quadrant = QUADRANTS(cornerAdjacentCell);
hasReturned = false;

if HAS_SAME_PARENT(quadrant, direction)
    % CASE 1: hasSameParent process
    cornerAdjacentCell = cornerAdjacentCell + DELTA_QUADRANT_HAS_SAME_PARENT(quadrant);
    while nodeState(cornerAdjacentCell) == -1
        cornerAdjacentCell = QTree_searchChildNode(cornerAdjacentCell, quadrant);
    end
else
    % CASE 2: isParentNodeAdjacent processs
    ancestorQuadrants = zeros(10,1); % WARNING: Need to change from 10 to MAX_TREE_DEPTH.
    ancestorQuadrants(1) = quadrant;
    indexQuadrant = 1;
    isNotParentNodeAdjacent = PARENT_NODE_ADJACENT_DIRECTION(quadrant, direction) == 0;

    if isNotParentNodeAdjacent
        % CASE 3: isNotParentNodeAdjacent processs
        isSameAncestor = true;
        while isSameAncestor
            indexQuadrant = indexQuadrant + 1;

            [cornerAdjacentCell, ancestorQuadrants(indexQuadrant)] = QTree_searchParentNode(cornerAdjacentCell, QUADRANTS);

            if cornerAdjacentCell == 1
                cornerAdjacentCell = 0;
                hasReturned = true;
                break
            end
            if HAS_SAME_PARENT(ancestorQuadrants(indexQuadrant), direction)
                cornerAdjacentCell = cornerAdjacentCell + DELTA_QUADRANT_HAS_SAME_PARENT(ancestorQuadrants(indexQuadrant));
                while nodeState(cornerAdjacentCell) == -1
                    cornerAdjacentCell = QTree_searchChildNode(cornerAdjacentCell, ancestorQuadrants(indexQuadrant));
                end
                hasReturned = true;
                break
            end
            if ancestorQuadrants(indexQuadrant) ~= ancestorQuadrants(indexQuadrant - 1) % Situation not changed. Need to more loop.
                isSameAncestor = false;
            end
        end
        if ~hasReturned && PARENT_NODE_ADJACENT_DIRECTION(ancestorQuadrants(indexQuadrant),direction) == 0
            hasReturned = true;
        end
    end

    if ~hasReturned
        [cornerAdjacentCell, ~] = QTree_searchParentNode(cornerAdjacentCell, QUADRANTS);

        if cornerAdjacentCell == 1
            cornerAdjacentCell = 0; % NOTE: No corner neighborhood node.
            hasReturned = true;
        end
    end

    if ~hasReturned
        [cornerAdjacentCell,hasPolygonInEdge] = QTree_findEdgeCell(cornerAdjacentCell,PARENT_NODE_ADJACENT_DIRECTION(ancestorQuadrants(indexQuadrant),direction),nodeState,QUADRANTS);

        if cornerAdjacentCell == 0
            % Nothing to do
        elseif isNotParentNodeAdjacent && hasPolygonInEdge % if cornerAdjacentCell is the leaf node
            cornerAdjacentCell = 0; % NOTE: No corner neighborhood node. 
        elseif nodeState(cornerAdjacentCell) == -1 % If neighborhood node is the internal node
            cornerAdjacentCell = QTree_searchChildNode(cornerAdjacentCell, SYMMETRIC_CORNER_DIRECTION(ancestorQuadrants(indexQuadrant)));
            while nodeState(cornerAdjacentCell) == -1
                cornerAdjacentCell = QTree_searchChildNode(cornerAdjacentCell, SYMMETRIC_CORNER_DIRECTION(direction));
            end
        end
    end
end

return