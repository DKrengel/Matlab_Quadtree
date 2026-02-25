function node = QTree_whichNode(nd, x_polygon, y_polygon, cellCenterX, cellCenterY)

% PURPOSE: 
%   Determine which child node of a quadtree cell a polygon belongs to.
% INPUT
%   nd                         = ID of the current (parent) node
%   x_polygon /y_polygon       = x / y -coordinate of the polygon
%   cellCenterX / cellCenterY  = x / y -coordinates of cell centers
% OUTPUT
%   node         = ID of the child node to which the polygon belongs
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

if cellCenterY(nd) <= y_polygon
    if x_polygon < cellCenterX(nd)
        node = 4 * nd - 2;
    else
        node = 4 * nd - 1;
    end
else
    if x_polygon < cellCenterX(nd)
        node = 4 * nd;
    else
        node = 4 * nd + 1;
    end
end

return
