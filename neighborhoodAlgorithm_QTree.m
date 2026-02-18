function [overlapPairs, overlapPairsCount] = neighborhoodAlgorithm_QTree(polygonVerticesCount, polygonX, polygonY, wallCount)

% Global variables for quadtree neighborhood algorithm.
persistent isInitialized MAX_TREE_DEPTH MAX_LEAF_NODES QUADRANTS fi_n la_n
persistent cellCenterX cellCenterY cellBottomX cellBottomY cellTopX cellTopY
persistent boxToNodeMap nodeState boxesCount
persistent boxBottomX boxTopX boxBottomY boxTopY boxCenterX boxCenterY xBoxMax yBoxMax
persistent minIndexArrayOfNeighborBoxPairs maxIndexArrayOfNeighborBoxPairs

% PURPOSE: 
% Neighborhood computation via a quadtree-approach
% INPUT: 
%  polygonVerticesCount = number of vertices per particle
%  polygonX, polygonY = vertex coordinates
%  wallCount = number of walls
% OUTPUT: 
%  overlapPairs = List of possible contacting polygons
%  overlapPairsCount = Number of possible contacting polygons
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

% For splitted bounding boxes
% persistent boundingBoxPolygonMapping

if isempty(isInitialized)
    [
        boxToNodeMap,nodeState,MAX_TREE_DEPTH,MAX_LEAF_NODES,QUADRANTS,cellCenterX,cellCenterY,...
        cellBottomX,cellBottomY,cellTopX,cellTopY,boxesCount,overlapPairs,overlapPairsCount,~,...
        boxBottomX, boxTopX, boxBottomY, boxTopY, boxCenterX, boxCenterY, xBoxMax, yBoxMax,...
        minIndexArrayOfNeighborBoxPairs, maxIndexArrayOfNeighborBoxPairs, fi_n, la_n...
    ] = QTree_initialize(polygonVerticesCount, polygonX, polygonY, wallCount);
    isInitialized = true;
    %return
end

%{
%Visualize the quadtree together with the system configuration
figure(3), clf, hold on
draw_tree(...
    MAX_LEAF_NODES, boxesCount, nodeState, boxCenterX, boxCenterY, ...
    cellCenterX, cellCenterY, min(boxBottomX), min(boxBottomY), ...
    max(boxTopX), max(boxTopY), MAX_TREE_DEPTH, boxTopX, boxTopY, ...
    boxBottomX, boxBottomY, fi_n, la_n, polygonVerticesCount, polygonX, polygonY,...
    cellBottomX, cellBottomY, cellTopX, cellTopY...
);
%}

[boxBottomX, boxTopX, boxBottomY, boxTopY, boxCenterX, boxCenterY] = ...
  QTree_updateBoundingBox(polygonVerticesCount, polygonX, polygonY, ...
  boxBottomX, boxTopX, boxBottomY, boxTopY);

hasChanges = false(boxesCount, 1);
isChanged  = false;

% Check whether boxes are overflow from each cell or not
isOverflow = boxCenterX < cellBottomX(boxToNodeMap) | cellTopX(boxToNodeMap) < boxCenterX | boxCenterY < cellBottomY(boxToNodeMap) | cellTopY(boxToNodeMap) < boxCenterY;
boxesToChange = find(isOverflow);
boxCount = length(boxesToChange);

for e = 1:boxCount
    baseBox = boxesToChange(e); % Base Box number
    hasChanges(baseBox) = true;
    baseNode = boxToNodeMap(baseBox); % Base node number
    % TODO: As it is, parallelization is not possible, so it is necessary to consider
    while nodeState(baseNode) == 0
        baseNode = (baseNode + 3 - QUADRANTS(baseNode))/4;
    end
    nodeState(baseNode) = 0; % Make the node "baseNode" empty

    % Make the bottom-up-process here
    switch QUADRANTS(baseNode)
        case 1
            otherNodes = [baseNode+1 baseNode+2 baseNode+3];
        case 2
            otherNodes = [baseNode-1 baseNode+1 baseNode+2];
        case 3
            otherNodes = [baseNode-2 baseNode-1 baseNode+1];
        case 4
            otherNodes = [baseNode-3 baseNode-2 baseNode-1];
    end

    leafOrInternalCount = 0;
    isNeededBottomUp = false;

    % Search Leaf or Internal node from the other 3 nodes
    for i = otherNodes
        isLeafOrInternal = QTree_isLeafOrInternalNode(i, nodeState);

        if isLeafOrInternal
            leafOrInternalCount = leafOrInternalCount + 1;
            nodeToBottomUp = i;
        end
    end

    % If only one node is a Leaf node among 3 nodes.
    if leafOrInternalCount == 1
        isLeaf = QTree_isLeafNode(nodeToBottomUp, nodeState);
        if isLeaf
            [newNode,~] = QTree_searchParentNode(nodeToBottomUp,QUADRANTS);
            boxToNodeMap(nodeState(nodeToBottomUp)) = newNode;
            nodeState(newNode) = nodeState(nodeToBottomUp);
            nodeState(nodeToBottomUp) = 0;
            isNeededBottomUp = true;
            hasChanges(nodeState(newNode)) = true;
        end
    end

    % If more bottom-up-process is needed
    while isNeededBottomUp
        switch QUADRANTS(newNode)
            case 1
                otherNodes = [newNode+1 newNode+2 newNode+3];
            case 2
                otherNodes = [newNode-1 newNode+1 newNode+2];
            case 3
                otherNodes = [newNode-2 newNode-1 newNode+1];
            case 4
                otherNodes = [newNode-3 newNode-2 newNode-1];
        end
        emptyCount = 0;
        % Search EMPTY node from the other 3 nodes
        for i = otherNodes
            isEmpty = QTree_isEmptyNode(i, nodeState);
            if isEmpty
                emptyCount = emptyCount + 1;
            end
        end
        % If all 3 nodes are EMPTY nodes
        if emptyCount == 3 % If 3 node around newNode is EMPTY node
            polygon_num = nodeState(newNode);
            nodeState(newNode) = 0;

            [newNode,~] = QTree_searchParentNode(newNode,QUADRANTS);

            boxToNodeMap(polygon_num) = newNode;
            nodeState(newNode) = polygon_num;
        else
            isNeededBottomUp = false; % Another bottom-up-process is not needed
        end
    end

    % Move node to the newly assigned node
    x = boxCenterX(baseBox);
    y = boxCenterY(baseBox);
    node = 1;
    for dep = 1:MAX_TREE_DEPTH % top-down-loop for decide the new place (node)
        node = QTree_whichNode(node, x, y, cellCenterX, cellCenterY);

        switch nodeState(node)
            case 0
                boxToNodeMap(baseBox) = node;
                nodeState(node) = baseBox;
                break
            case -1
                continue
            otherwise
                nodeNext = node;
                boxConflict = find(boxToNodeMap == node);
                hasChanges(boxConflict) = true;
                while true
                    nodeState(nodeNext) = -1;
                    nodeNextConflict = QTree_whichNode(nodeNext, ...
                        boxCenterX(boxConflict), boxCenterY(boxConflict), cellCenterX, cellCenterY);

                    nodeNextOriginal = QTree_whichNode(nodeNext, ...
                        x, y, cellCenterX, cellCenterY);

                    if nodeNextConflict == nodeNextOriginal
                        nodeNext = nodeNextConflict;
                    else
                        boxToNodeMap(baseBox)       = nodeNextOriginal;
                        boxToNodeMap(boxConflict)   = nodeNextConflict;
                        nodeState(nodeNextOriginal) = baseBox;
                        nodeState(nodeNextConflict) = boxConflict;
                        isChanged = true;
                        break
                    end
                end
                if isChanged
                    break
                end
        end
    end
end

needsRevisiting = find(hasChanges); % Column vector
if ~isempty(needsRevisiting)
    recheckNodes  = boxToNodeMap(needsRevisiting)';
    maxPairsCount = 3 * boxesCount;
    sumBoxes1 = zeros(maxPairsCount, 1);
    sumBoxes2 = zeros(maxPairsCount, 1);
    index = 0;
    for cell = recheckNodes
        [boxes1, boxes2, boxesLength] = QTree_searchNeighborBoxesPairDiffs(cell, nodeState, QUADRANTS, maxPairsCount,...
        cellCenterX, cellCenterY, cellBottomX, cellBottomY, cellTopX, cellTopY, xBoxMax, yBoxMax, fi_n, la_n);

        sumBoxes1((index + 1):(index + boxesLength)) = boxes1;
        sumBoxes2((index + 1):(index + boxesLength)) = boxes2;
        index = index + boxesLength;
    end
    isUnchangedMin = ~ismember(minIndexArrayOfNeighborBoxPairs, needsRevisiting);
    isUnchangedMax = ~ismember(maxIndexArrayOfNeighborBoxPairs, needsRevisiting);
    isUnchangedBoth = isUnchangedMin & isUnchangedMax;
    % Update min/maxIndexArrayOfNeighborBoxPairs
    minIndexArrayOfNeighborBoxPairs = [ minIndexArrayOfNeighborBoxPairs(isUnchangedBoth); sumBoxes1(1:index) ];
    maxIndexArrayOfNeighborBoxPairs = [ maxIndexArrayOfNeighborBoxPairs(isUnchangedBoth); sumBoxes2(1:index) ];
end

if nnz(minIndexArrayOfNeighborBoxPairs) > 0
    % % For splitted bounding boxes
    % isUnique = find(boundingBoxPolygonMapping(minIndexArrayOfNeighborBoxPairs) - boundingBoxPolygonMapping(b2));
    % minIndexArrayOfNeighborBoxPairs = minIndexArrayOfNeighborBoxPairs(isUnique);
    % maxIndexArrayOfNeighborBoxPairs = maxIndexArrayOfNeighborBoxPairs(isUnique);

    isOverlapped = boxTopY(minIndexArrayOfNeighborBoxPairs) > boxBottomY(maxIndexArrayOfNeighborBoxPairs) & ...
                   boxBottomY(minIndexArrayOfNeighborBoxPairs) < boxTopY(maxIndexArrayOfNeighborBoxPairs) & ...
                   boxTopX(minIndexArrayOfNeighborBoxPairs) > boxBottomX(maxIndexArrayOfNeighborBoxPairs) & ...
                   boxBottomX(minIndexArrayOfNeighborBoxPairs) < boxTopX(maxIndexArrayOfNeighborBoxPairs);
    col1 = minIndexArrayOfNeighborBoxPairs(isOverlapped);
    col2 = maxIndexArrayOfNeighborBoxPairs(isOverlapped);

    % % For splitted bounding boxes
    % col1 = boundingBoxPolygonMapping(col1);
    % col2 = boundingBoxPolygonMapping(col2);

    overlapPairs = [col1, col2]';
    if isempty(overlapPairs)
        overlapPairs = zeros(2,1);
    end
else
    overlapPairs = zeros(2,1);
end

overlapPairsCount = size(overlapPairs,2);

return