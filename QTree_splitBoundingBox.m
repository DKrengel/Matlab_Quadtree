function [topX, bottomX, topY, bottomY, boundingBoxCount, isXAxisAligned] =...
           QTree_splitBoundingBox(x, y, splitSize)
% PURPOSE: 
% split bounding boxes for large particles (i.e. walls)
% INPUT: 
%  x, y      = vertex coordinates
%  splitSize = split size
% OUTPUT: 
%  topX, bottomX, topY, bottomY = Top / Bottom x / y coordinates of bounding box
%  boundingBoxCount  = total number of bounding boxes
%  isXAxisAligned    = whether split bounding boxes are aligned with the x-axis or not
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

MAX_SPLIT_COUNT = 100;

topX    = NaN(1,MAX_SPLIT_COUNT);
bottomX = NaN(1,MAX_SPLIT_COUNT);
topY    = NaN(1,MAX_SPLIT_COUNT);
bottomY = NaN(1,MAX_SPLIT_COUNT);

xMax = max(x);
xMin = min(x);
yMax = max(y);
yMin = min(y);
xLength = xMax - xMin;
yLength = yMax - yMin;

if xLength > yLength
    isXAxisAligned = true;
else
    isXAxisAligned = false;
end

[a, b, c] = getEdgeEquation(x, y, isXAxisAligned);

pastXIntersection = [];
pastYIntersection = [];
isBreak = false;
if isXAxisAligned
    goal = xMin;
else
    goal = yMin;
end
for iSplit = 1:MAX_SPLIT_COUNT
    if isXAxisAligned
        start = goal;
        goal = start + splitSize;
        yCandidate = - (a .* goal + c) ./ b;
        candidateIntersectionNum = length(yCandidate);
        xCandidate = zeros(1,candidateIntersectionNum) + goal;
        isInside   = pointsOnEdge(xCandidate,yCandidate,x,y);
        if nnz(isInside) > 0
            xIntersection = xCandidate(isInside);
            yIntersection = yCandidate(isInside);
            n = (start <= x) & (x < goal);
            xCornerCandidate = [xIntersection pastXIntersection x(n)];
            yCornerCandidate = [yIntersection pastYIntersection y(n)];
        else
            n = (start <= x) & (x <= xMax);
            xCornerCandidate = [pastXIntersection x(n)];
            yCornerCandidate = [pastYIntersection y(n)];
            isBreak = true;
        end
    else
        start = goal;
        goal = start + splitSize;
        xCandidate = - (b .* goal + c) ./ a;
        candidateIntersectionNum = length(xCandidate);
        yCandidate = zeros(1,candidateIntersectionNum) + goal;
        isInside   = pointsOnEdge(xCandidate,yCandidate,x,y);
        if nnz(isInside) > 0
            xIntersection = xCandidate(isInside);
            yIntersection = yCandidate(isInside);
            n = (start <= y) & (y < goal);
            xCornerCandidate = [xIntersection pastXIntersection x(n)];
            yCornerCandidate = [yIntersection pastYIntersection y(n)];
        else
            n = (start <= y) & (y <= yMax);
            xCornerCandidate = [pastXIntersection x(n)];
            yCornerCandidate = [pastYIntersection y(n)];
            isBreak = true;
        end
    end
    topX(iSplit)    = max(xCornerCandidate);
    bottomX(iSplit) = min(xCornerCandidate);
    topY(iSplit)    = max(yCornerCandidate);
    bottomY(iSplit) = min(yCornerCandidate);
    if isBreak
        break
    end
    pastXIntersection = xIntersection;
    pastYIntersection = yIntersection;
end

boundingBoxCount = iSplit;

end

%{
    Function of extended inpolygon.
    ref: https://jp.mathworks.com/matlabcentral/answers/285565-plot
%}
function pointsOnEdge = pointsOnEdge(x, y, xv, yv)
%PURPOSE: get the equation for each edge of the polyon
%         ax + by + c = 0
ERROR = 1.005; %error tolerance

%  @return boundingBoxCount {double} Number of bounding boxes
xvExpanded = ERROR * xv + (1 - ERROR) * mean(xv(1:end-1));
yvExpanded = ERROR * yv + (1 - ERROR) * mean(yv(1:end-1));

pointsOnEdge = inpolygon(x, y, xvExpanded, yvExpanded);

%{
    Function of get each edge equation of polygon.
    Equation: ax + by + c = 0

    @param x {double(1, #corner)} X coordinates of wall corners
    @param y {double(1, #corner)} Y coordinates of wall corners

    @returns a {double(1, #edge)} X coefficients
    @returns b {double(1, #edge)} Y coefficients
    @returns c {double(1, #edge)} Intercept
%}

end



function [a, b, c] = getEdgeEquation(x, y, isXAxisAligned)
% PURPOSE:

 x1 = x;
 x2 = [x(2:end) x(1)];
 y1 = y;
 y2 = [y(2:end) y(1)];
 a = (y1 - y2);
 b = - (x1 - x2);
 c = (y1 - y2) .* (-x1) + y1 .* (x1 - x2);

 % Deletion of cases where straight lines do not intersect.
 if isXAxisAligned
    n = b == 0;
 else
    n = a == 0;
 end
 a(:,n) = [];
 b(:,n) = [];
 c(:,n) = [];
end
