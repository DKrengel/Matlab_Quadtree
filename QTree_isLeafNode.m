function isLeaf = QTree_isLeafNode(nodeToBottomUp, nodeState)

% PURPOSE: 
% Decide whether node i is leaf node
% INPUT: 
%  nodeToBottomUp = 
%  nodeState
% OUTPUT: 
%  isLeaf
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
