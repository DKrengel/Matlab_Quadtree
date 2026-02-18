function [boxBottomX, boxTopX, boxBottomY, boxTopY, boxCenterX, boxCenterY] = ...
           QTree_updateBoundingBox(polygonVerticesCount, polygonX, polygonY, ...
           boxBottomX, boxTopX, boxBottomY, boxTopY)

% PURPOSE: 
% INPUT: 
%  polygonVerticesCount = number of vertices per particle
%  polygonX, polygonY = vertex coordinates
%  boxBottomX, boxTopX, boxBottomY, boxTopY
% OUTPUT: 
%  boxBottomX, boxTopX = X coordinates of the bottom / top of the bounding box
%  boxBottomY, boxTopY = Y coordinates of the bottom / top of the bounding box
%  boxCenterX, boxCenterY = X / Y coordinates of the center of the bounding box
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

% Check if it is necessary to calculate every bounding box of the wall

boundingBoxCount = length(boxBottomX);
polygonCount = length(polygonVerticesCount);
bRange = 1:boundingBoxCount;
pRange = 1:polygonCount;

% return {double(1, #polygons)}
[polygonBottomX,polygonTopX] = update_boundbox(polygonVerticesCount(pRange),polygonX(:, pRange));
[polygonBottomY,polygonTopY] = update_boundbox(polygonVerticesCount(pRange),polygonY(:, pRange));

boxBottomX(bRange) = polygonBottomX;
boxTopX(bRange) = polygonTopX;
boxBottomY(bRange) = polygonBottomY;
boxTopY(bRange) = polygonTopY;

boxCenterX = (boxBottomX + boxTopX) ./ 2;
boxCenterY = (boxBottomY + boxTopY) ./ 2;


return