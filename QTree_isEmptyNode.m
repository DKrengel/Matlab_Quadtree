function isEmpty = QTree_isEmptyNode(i, nodeState)

% PURPOSE: 
% Decide whether node i is empty
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

isEmpty = false;

if nodeState(i) == 0
    isEmpty = true;
end

return