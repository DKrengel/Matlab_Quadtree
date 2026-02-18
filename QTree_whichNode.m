function node = QTree_whichNode(nd, x_polygon, y_polygon, cellCenterX, cellCenterY)

% PURPOSE: 
% INPUT: 
%  nd = 
%  x_polygon  y_polygon
% cellCenterX / cellCenterY
% OUTPUT: 
%  node
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
