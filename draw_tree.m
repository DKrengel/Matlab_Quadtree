function draw_tree(...
    MAX_LEAF_NODES, boxesCount, nodeState, boxCenterX, boxCenterY,...
    cellCenterX, cellCenterY, X_SYSTEM_MIN, Y_SYSTEM_MIN,...
    X_SYSTEM_MAX, Y_SYSTEM_MAX, MAX_TREE_DEPTH, boxTopX, boxTopY,...
    boxBottomX, boxBottomY, fi_n, la_n, polygonVerticesCount, polygonX, polygonY,...
    cellBottomX, cellBottomY, cellTopX, cellTopY...
)

%Purpose:
% draw the 3D tre structure over a 2D view of the system

DEF_MARKER = 20; % (SYSTEM_SIZE=16):30
MARKER_PARAMETER = 0.8;
DEF_LINE_WIDTH = 1;
HEIGHT = 0; % NOTE: How high from the bottom the tree should be displayed when plotting

plots = zeros(boxesCount, 5);

% Color
cFace = [0.8 0.7 0.8];
cEdge = [0.5 0.2 0.55];
cPolygon = [160.0/255.0 160.0/255.0 160.0/255.0];

% Plot funnel
n=length(polygonVerticesCount);

for i = 1:n
    nVertex = polygonVerticesCount(i);
    plots(i, 1) = patch(polygonX(1:nVertex, i), polygonY(1:nVertex, i), zeros(nVertex, 1) + HEIGHT, cPolygon);
end

for iPlot = 1:boxesCount
    tx = boxTopX(iPlot);
    bx = boxBottomX(iPlot);
    ty = boxTopY(iPlot);
    by = boxBottomY(iPlot);
    plots(i, 2) = patch([bx tx tx bx], [by by ty ty], zeros(1, 4) + HEIGHT, cFace, 'LineWidth', 1, 'EdgeColor', cEdge);
    cx = boxCenterX(iPlot);
    cy = boxCenterY(iPlot);
    plot3(cx, cy, 0 + HEIGHT, 'k*', 'MarkerSize', 4)
end

alpha(0.5)

% ===plot tree===
intl_feature = zeros(MAX_LEAF_NODES-boxesCount,1);
leaf_feature = zeros(boxesCount,1);
j=1;
for i = 1:MAX_LEAF_NODES % node-number
    num = nodeState(i); % polygon-number
    switch num
        case 0
            continue
        case -1
            intl_feature(j) = i;
            j = j+1;
        otherwise
            leaf_feature(num) = i;
    end
end

plot3(0,0,0,'MarkerSize', .1) % Define 3d plot area before hold on.

tr_plot_y = zeros(boxesCount,1);
plots_tree = zeros(MAX_LEAF_NODES, 4);

diff_dep_node = -0.1; % (SYSTEM_SIZE=16):-1.5
for i = 1:MAX_TREE_DEPTH
    j = fi_n(i):la_n(i);
    tr_plot_y(j) = diff_dep_node*i;
end
height_node_z = 1; % (SYSTEM_SIZE=16):21 % MEMO: Shortened the distance between tree and system
tr_plot_y = tr_plot_y + height_node_z;

% Plot line from node to node
depth = 1;
for i = 2:MAX_LEAF_NODES
    if nodeState(i)~=0
        yy=[tr_plot_y(i) tr_plot_y(floor((i+2)/4))];
        if nodeState(i)==-1
            xx=[cellCenterY(i) cellCenterY(floor((i+2)/4))];
            zz=[cellCenterX(i) cellCenterX(floor((i+2)/4))];
            plots_tree(i,1) = plot3(zz, xx, yy, 'k-', 'LineWidth', DEF_LINE_WIDTH*(1/depth));
        else
            ii = nodeState(i);
            xx=[boxCenterX(ii) cellCenterX(floor((i+2)/4))];
            zz=[boxCenterY(ii) cellCenterY(floor((i+2)/4))];
            plots_tree(i,1) = plot3(xx, zz, yy, 'k-');
        end
    end
    if la_n(depth) == i
        depth = depth + 1;
    end
end

% Plot Internal node
% MEMO: Using for-loop to change 'MarkerSize' by tree height.
internal_nodes = intl_feature(intl_feature~=0)';
depth = 1;
for i = internal_nodes
    plots_tree(i,2) = plot3(cellCenterX(i), cellCenterY(i), tr_plot_y(i), 'ko', 'MarkerSize', DEF_MARKER*MARKER_PARAMETER^(depth-1), 'MarkerFaceColor','white', 'LineWidth',1);
    if la_n(depth) == i
        depth = depth + 1;
    end
end

% MEMO: Faster than for-loop but cannot change 'MarkerSize' by tree height.
% i = intl_feature(intl_feature~=0);
% plots_tree(i,2) = plot3(cellCenterX(i), cellCenterY(i), tr_plot_y(i), 'ko', 'MarkerSize', DEF_MARKER, 'MarkerFaceColor','white', 'LineWidth',1);

% Plot leaf node
dl = 0.005; % (SYSTEM_SIZE=16):.3
for i = 1:length(leaf_feature)
    ii = leaf_feature(i);
    % disp(['polygon: ' num2str(i) ', node: ' num2str(ii)]) % node-number
    zcen = tr_plot_y(ii);
    nd_ii = nodeState(ii);
    xcen = boxCenterX(nd_ii);
    ycen = boxCenterY(nd_ii);
    plots_tree(i,3) = plot3([xcen-dl xcen+dl xcen+dl xcen-dl xcen-dl], [ycen ycen ycen ycen ycen], [zcen+dl zcen+dl zcen-dl zcen-dl zcen+dl],'b');
end

% Plot Line between leaf-node and polygon
% plots_line = zeros(length(leaf_feature),1);
% for i = 1:length(leaf_feature)
%     ii = leaf_feature(i);
%     plots_line(i,1) = plot3([boxCenterX(i) boxCenterX(i)], [boxCenterY(i) boxCenterY(i)], [HEIGHT tr_plot_y(ii)], 'k-', 'LineWidth', DEF_LINE_WIDTH/10); % [polygon tree-node]
% end

% NOTE: Plot each cell
plot3([X_SYSTEM_MIN X_SYSTEM_MIN X_SYSTEM_MAX X_SYSTEM_MAX X_SYSTEM_MIN],[Y_SYSTEM_MAX Y_SYSTEM_MIN Y_SYSTEM_MIN Y_SYSTEM_MAX Y_SYSTEM_MAX],[0 0 0 0 0]+HEIGHT, 'k-')
plotCells = zeros(length(intl_feature(intl_feature~=0)), 2);
for i = 1:length(intl_feature(intl_feature~=0))
    ii = intl_feature(i);
    plotCells(i,1) = plot3([cellBottomX(ii) cellTopX(ii)], [cellCenterY(ii) cellCenterY(ii)], [0 0]+HEIGHT, 'k-');
    plotCells(i,2) = plot3([cellCenterX(ii) cellCenterX(ii)], [cellBottomY(ii) cellTopY(ii)], [0 0]+HEIGHT, 'k-');
end

hold off
xlim([X_SYSTEM_MIN X_SYSTEM_MAX])
ylim([Y_SYSTEM_MIN Y_SYSTEM_MAX])
zlim([HEIGHT max(tr_plot_y)])
axis image
view([-15 15])
zticks([])

% MEMO: Turn off z-axis
Ax = gca;
Ax.ZAxis.Visible = 'off';
Ax.ZGrid = 'off';
Ax.Color = 'none';
