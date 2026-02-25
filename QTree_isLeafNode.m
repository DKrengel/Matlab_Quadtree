function isLeaf = QTree_isLeafNode(nodeToBottomUp, nodeState)

% PURPOSE: 
% Decide whether node i is leaf node
% INPUT: 
%   nodeToBottomUp = node index to be checked
%   nodeState      = node state array (empty/internal/leaf-with-box-id)
% OUTPUT
%   isLeaf         = logical value indicating whether the node is a leaf node
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

isLeaf = false;

if nodeState(nodeToBottomUp) > 0
    isLeaf = true;
end

return
