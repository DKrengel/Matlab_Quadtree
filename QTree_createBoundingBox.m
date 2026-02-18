function [bottomX, topX, bottomY, topY, centerX, centerY, ...
          polygonMapping, pMin, xBoxMax, yBoxMax] = QTree_createBoundingBox(...
          polygonVerticesCount, polygonX, polygonY, wallCount)

% PURPOSE: 
% create the particle bounding boxes suitable for the quadtree
% INPUT: 
%  polygonVerticesCount = number of vertices per particle
%  polygonX, polygonY = vertex coordinates
%  wallCount = number of walls
% OUTPUT: 
% bottomX, topX = X coordinates of the bottom / top of the bounding box
% bottomY, topY = Y coordinates of the bottom / top of the bounding box
% centerX, centerY = X / Y coordinates of the center of the bounding box
% polygonMapping   = Polygon number of each bounding box
% pMin             = Minimum particle size (width or height)
% xBoxMax, yBoxMax = ??
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori


SIZE_SCALING_FACTOR = 100;

polygonCount = length(polygonVerticesCount);
maxBoundingBoxCount = polygonCount * 5;

% Definition of output values
bottomX = zeros(1, maxBoundingBoxCount);
topX    = zeros(1, maxBoundingBoxCount);
bottomY = zeros(1, maxBoundingBoxCount);
topY    = zeros(1, maxBoundingBoxCount);
polygonMapping = zeros(maxBoundingBoxCount, 1);

[polygonBottomX,polygonTopX] = update_boundbox(polygonVerticesCount,polygonX);
[polygonBottomY,polygonTopY] = update_boundbox(polygonVerticesCount,polygonY);

w = polygonTopX - polygonBottomX;
h = polygonTopY - polygonBottomY;
start = wallCount + 1;

xBoxMax = max(w(1:end)); % including wall. if do not want to include, change to w(start:end)
yBoxMax = max(h(1:end)); % including wall. if do not want to include, change to h(start:end)

pMin = min([w(start:end) h(start:end)]);
splitSize = pMin * SIZE_SCALING_FACTOR;
boxSize = max([w; h]); % (1, #polygons)

lastIndex = 1;
for iPolygon = 1:polygonCount
    if boxSize(iPolygon) < splitSize
        bottomX(lastIndex) = polygonBottomX(iPolygon);
        topX(lastIndex)    = polygonTopX(iPolygon);
        bottomY(lastIndex) = polygonBottomY(iPolygon);
        topY(lastIndex)    = polygonTopY(iPolygon);
        polygonMapping(lastIndex) = iPolygon;
        lastIndex = lastIndex + 1;
    else
        vertex = polygonVerticesCount(iPolygon);
        x = polygonX(1:vertex, iPolygon)';
        y = polygonY(1:vertex, iPolygon)';
        [tx, bx, ty, by, n, ~] = QTree_splitBoundingBox(x, y, splitSize); % 'n' must be 2 or more
        iRange = 1:n;
        oRange = lastIndex:lastIndex+n-1;
        bottomX(oRange) = bx(iRange);
        topX(oRange)    = tx(iRange);
        bottomY(oRange) = by(iRange);
        topY(oRange)    = ty(iRange);
        polygonMapping(oRange) = iPolygon;
        lastIndex = lastIndex+n;
    end
end

splitedBoundingBoxCount = lastIndex - 1;
range = 1:splitedBoundingBoxCount;

bottomX = bottomX(range);
topX    = topX(range);
bottomY = bottomY(range);
topY    = topY(range);
polygonMapping = polygonMapping(range);

centerX = (bottomX + topX) ./ 2;
centerY = (bottomY + topY) ./ 2;

return