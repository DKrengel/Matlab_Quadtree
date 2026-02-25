function isEmpty = QTree_isEmptyNode(i, nodeState)

% PURPOSE: 
% Decide whether node i is empty
% INPUT: 
%   i         = node index to be checked
%   nodeState = node state array (empty/internal/leaf-with-box-id)
% OUTPUT: 
%   isEmpty   = logical value indicating whether the node is empty
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