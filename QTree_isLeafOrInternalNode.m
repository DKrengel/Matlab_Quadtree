function isLeafOrInternal = QTree_isLeafOrInternalNode(i, nodeState)

% PURPOSE: 
% Decide whether node i is leaf or internal / branch node
% INPUT: 
%  i = 
%  nodeState
% OUTPUT: 
%  isLeafOrInternal
% CAVEAT: 
% PERFORMANCE CONSIDERATION: 
% TODO: 
% LITERATURE: 
% REVISION HISTORY: Y.Watanabe & K.Kandori

isLeafOrInternal = true;

% Definition of nodeState is
% leaf node: larger than 0
% Empty node: 0
% Internal node: -1
if nodeState(i) == 0 % If empty
    isLeafOrInternal = false;
end

return