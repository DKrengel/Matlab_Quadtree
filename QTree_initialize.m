function [...
    boxToNodeMap,nodeState,MAX_TREE_DEPTH,MAX_LEAF_NODES,QUADRANTS,cellCenterX,cellCenterY,...
    cellBottomX,cellBottomY,cellTopX,cellTopY,boxesCount,overlapPairs,overlapPairsCount,boundingBoxPolygonMapping,...
    boxBottomX, boxTopX, boxBottomY, boxTopY, boxCenterX, boxCenterY, xBoxMax, yBoxMax,...
    minIndexArrayOfNeighborBoxPairs, maxIndexArrayOfNeighborBoxPairs, fi_n, la_n...
] = QTree_initialize(polygonVerticesCount, polygonX, polygonY, wallCount)

% PURPOSE: 
% Initialize the quadtree at first use
% INPUT: 
%  polygonVerticesCount = number of vertices per particle
%  polygonX, polygonY = vertex coordinates
%  wallCount = number of walls
% OUTPUT: 
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

%create the bounding boxes for the particles; split the bounding box where
% necessary into multiples (i.e. for long walls)
[boxBottomX, boxTopX, boxBottomY, boxTopY, boxCenterX, boxCenterY,...
 boundingBoxPolygonMapping, minimumParticleSize, xBoxMax, yBoxMax] = ...
 QTree_createBoundingBox(polygonVerticesCount, polygonX, polygonY, wallCount);

boxesCount = length(boundingBoxPolygonMapping);
X_SYSTEM_MIN = min(boxBottomX);
X_SYSTEM_MAX = max(boxTopX);
Y_SYSTEM_MIN = min(boxBottomY);
Y_SYSTEM_MAX = max(boxTopY);

X_SYSTEM_SIZE = X_SYSTEM_MAX - X_SYSTEM_MIN;
Y_SYSTEM_SIZE = Y_SYSTEM_MAX - Y_SYSTEM_MIN;

divNo = max([X_SYSTEM_SIZE Y_SYSTEM_SIZE]) / minimumParticleSize;
MAX_TREE_DEPTH = round(log2(divNo)) + 2; 
boxToNodeMap = zeros(boxesCount, 1);

% Arrays of first or last node number
fi_n = ones(MAX_TREE_DEPTH, 1);
la_n = ones(MAX_TREE_DEPTH, 1);
sum_n = 1; % Variable to compute fi_n and la_n
for i = 2:MAX_TREE_DEPTH
    fi_n(i) = sum_n + 4^(i-2);
    la_n(i) = fi_n(i) + 4^(i-1)-1;
    sum_n = sum_n + 4^(i-2);
end
intl_max       = la_n(MAX_TREE_DEPTH-1);    % Last internal-node number
MAX_LEAF_NODES = la_n(MAX_TREE_DEPTH);      % Last leaf-node number
nodeState      = zeros(MAX_LEAF_NODES, 1);  % Empty-node: 0, Internal-node: -1, Leaf-node: polygon number

disp(['Max tree depth: ' num2str(MAX_TREE_DEPTH) ', Max leaf nodes: ' num2str(MAX_LEAF_NODES)])

cellCenterX = zeros(1, MAX_LEAF_NODES);
cellCenterY = zeros(1, MAX_LEAF_NODES);
cellBottomX = zeros(1, MAX_LEAF_NODES);
cellBottomY = zeros(1, MAX_LEAF_NODES);
cellTopX = zeros(1, MAX_LEAF_NODES);
cellTopY = zeros(1, MAX_LEAF_NODES);
x_delt_div = zeros(MAX_TREE_DEPTH, 1);
y_delt_div = zeros(MAX_TREE_DEPTH, 1);

for i = 1:MAX_TREE_DEPTH
    x_delt_div(i) = X_SYSTEM_SIZE/(2^i);
    y_delt_div(i) = Y_SYSTEM_SIZE/(2^i);
end
cellCenterX(1) = (X_SYSTEM_MIN + X_SYSTEM_MAX) / 2;
cellCenterY(1) = (Y_SYSTEM_MIN + Y_SYSTEM_MAX) / 2;
cellBottomX(1) = cellCenterX(1) - x_delt_div(1);
cellBottomY(1) = cellCenterY(1) - y_delt_div(1);
cellTopX(1) = cellCenterX(1) + x_delt_div(1);
cellTopY(1) = cellCenterY(1) + y_delt_div(1);

% Boundary coordinates of subdivisions
for i = 1:intl_max
    % Delta
    k = 0;
    tmp = i;
    while tmp > 4^k
        tmp = tmp - 4^k;
        k = k + 1;
    end
    x_delt = x_delt_div(k+2);
    y_delt = y_delt_div(k+2);
    x_mi_delt = [-x_delt x_delt -x_delt x_delt];
    y_mi_delt = [y_delt y_delt -y_delt -y_delt];
    % Middle
    cellCenterX(4*i-2:4*i+1) = cellCenterX(i) + x_mi_delt;
    cellCenterY(4*i-2:4*i+1) = cellCenterY(i) + y_mi_delt;
    
    % Lower
    cellBottomX(4*i-2:4*i+1) = cellCenterX(4*i-2:4*i+1) - x_delt;
    cellBottomY(4*i-2:4*i+1) = cellCenterY(4*i-2:4*i+1) - y_delt;
    
    % Upper
    cellTopX(4*i-2:4*i+1) = cellCenterX(4*i-2:4*i+1) + x_delt;
    cellTopY(4*i-2:4*i+1) = cellCenterY(4*i-2:4*i+1) + y_delt;
end

cnt_clas_part = 0;
depth  = 1;
node   = 1;
dep_od = true;
ma_len_ch = 0;
nd_gp_od  = (1:boxesCount)';

while cnt_clas_part < boxesCount
    % change the depth
    if node == fi_n(depth+1)
        depth  = depth + 1;
        dep_od = ~dep_od;
        if dep_od
            nd_gp_ev = zeros(ma_len_ch, 4^depth);
        else
            nd_gp_od = zeros(ma_len_ch, 4^depth);
        end
        ma_len_ch = 0;
    end
    if dep_od
        nd_clus = nd_gp_od(:, node-fi_n(depth)+1);
    else
        nd_clus = nd_gp_ev(:, node-fi_n(depth)+1);
    end
    nd_clus = nd_clus(nd_clus ~= 0);
    if length(nd_clus) > 1
        xn = cellCenterX(node);
        yn = cellCenterY(node);
        nw = nd_clus(boxCenterX(nd_clus)<=xn & boxCenterY(nd_clus)>yn);
        ne = nd_clus(boxCenterX(nd_clus)>xn & boxCenterY(nd_clus)>yn);
        sw = nd_clus(boxCenterX(nd_clus)<=xn & boxCenterY(nd_clus)<=yn);
        se = nd_clus(boxCenterX(nd_clus)>xn & boxCenterY(nd_clus)<=yn);
        nodeState(node) = -1;
        if length(nw) == 1
            boxToNodeMap(nw) = 4*node - 2;
            nodeState(4*node-2) = nw;
            cnt_clas_part = cnt_clas_part+1;
        end
        if length(ne) == 1
            boxToNodeMap(ne) = 4*node - 1;
            nodeState(4*node-1) = ne;
            cnt_clas_part = cnt_clas_part+1;
        end
        if length(sw) == 1
            boxToNodeMap(sw) = 4*node;
            nodeState(4*node) = sw;
            cnt_clas_part = cnt_clas_part+1;
        end
        if length(se) == 1
            boxToNodeMap(se) = 4*node + 1;
            nodeState(4*node+1) = se;
            cnt_clas_part = cnt_clas_part+1;
        end
        max_tmp = [length(nw); length(ne); length(sw); length(se)];
        if max(max_tmp) > ma_len_ch
            ma_len_ch = max(max_tmp);
        end
    end
    if dep_od
        nd_gp_ev(1:length(nw), 4*node-2-fi_n(depth+1)+1) = nw;
        nd_gp_ev(1:length(ne), 4*node-1-fi_n(depth+1)+1) = ne;
        nd_gp_ev(1:length(sw), 4*node-fi_n(depth+1)+1) = sw;
        nd_gp_ev(1:length(se), 4*node+1-fi_n(depth+1)+1) = se;
    else
        nd_gp_od(1:length(nw), 4*node-2-fi_n(depth+1)+1) = nw;
        nd_gp_od(1:length(ne), 4*node-1-fi_n(depth+1)+1) = ne;
        nd_gp_od(1:length(sw), 4*node-fi_n(depth+1)+1) = sw;
        nd_gp_od(1:length(se), 4*node+1-fi_n(depth+1)+1) = se;
    end
    nw=0; ne=0; sw=0; se=0;
    node = node + 1;
end

mtx = repmat([1 2 3 4]', 1, intl_max);
vec = reshape(mtx, [], 1);
QUADRANTS = [1;vec];

% Find neighboring pairs of the bounding boxes
maxPairsCount = 30 * boxesCount;
sumBoxes1 = zeros(maxPairsCount, 1);
sumBoxes2 = zeros(maxPairsCount, 1);
index = 0;
for iBox = 1:boxesCount
    iCell = boxToNodeMap(iBox);
    [boxes1, boxes2, bLen] = QTree_searchNeighborBoxes(iCell, nodeState, QUADRANTS, maxPairsCount,...
    cellCenterX, cellCenterY, cellBottomX, cellBottomY, cellTopX, cellTopY, xBoxMax, yBoxMax);

    sumBoxes1((index + 1):(index + bLen)) = boxes1;
    sumBoxes2((index + 1):(index + bLen)) = boxes2;
    index = index + bLen;
end

minIndexArrayOfNeighborBoxPairs = sumBoxes1(1:index);
maxIndexArrayOfNeighborBoxPairs = sumBoxes2(1:index);

if nnz(minIndexArrayOfNeighborBoxPairs) > 0
    % % For splitted bounding boxes
    % isUnique = find(boundingBoxPolygonMapping(minIndexArrayOfNeighborBoxPairs) - boundingBoxPolygonMapping(maxIndexArrayOfNeighborBoxPairs));
    % minIndexArrayOfNeighborBoxPairs = minIndexArrayOfNeighborBoxPairs(isUnique);
    % maxIndexArrayOfNeighborBoxPairs = maxIndexArrayOfNeighborBoxPairs(isUnique);

    isOverlapped = boxTopY(minIndexArrayOfNeighborBoxPairs) > boxBottomY(maxIndexArrayOfNeighborBoxPairs) & boxBottomY(minIndexArrayOfNeighborBoxPairs) < boxTopY(maxIndexArrayOfNeighborBoxPairs) & boxTopX(minIndexArrayOfNeighborBoxPairs) > boxBottomX(maxIndexArrayOfNeighborBoxPairs) & boxBottomX(minIndexArrayOfNeighborBoxPairs) < boxTopX(maxIndexArrayOfNeighborBoxPairs);
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