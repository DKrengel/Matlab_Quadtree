function [nwr,fp,sumArea] = weighted_normal(cut1x, cut1y, area1, cut2x, cut2y, area2, c1x, c1y, xCutCtr1,yCutCtr1, xCutCtr2,yCutCtr2)

tol   = 1.0e-10;
tolsq = tol^2;

if ( area1 < 0.0 || area2 < 0.0 )
   error('One of the area has negative sign')
end

% Find common point.
if ( (cut1x(1)-cut2x(1))^2+(cut1y(1)-cut2y(1))^2 < tolsq )
   id1 = [1 2];
   id2 = [1 2];
elseif ( (cut1x(1)-cut2x(2))^2+(cut1y(1)-cut2y(2))^2 < tolsq )
   id1 = [1 2];
   id2 = [2 1];
elseif ( (cut1x(2)-cut2x(2))^2+(cut1y(2)-cut2y(2))^2 < tolsq )
   id1 = [2 1];
   id2 = [2 1];
elseif ( (cut1x(2)-cut2x(1))^2+(cut1y(2)-cut2y(1))^2 < tolsq )
   id1 = [2 1];
   id2 = [1 2];
else
   error('Could not find common point.')
end
cut1x = cut1x(id1);
cut1y = cut1y(id1);
cut2x = cut2x(id2);
cut2y = cut2y(id2);

sumArea = area1 + area2;

% Weighted tangential direction.
r1 = [cut1x(2)-cut1x(1);cut1y(2)-cut1y(1)];
r2 = [cut2x(2)-cut2x(1);cut2y(2)-cut2y(1)];

wr = (r1*area1 + r2*area2) / sumArea;

% Check the direction of wr to obtain 
% the correct direction of the repulsive force.
vecar = [cut1x(1)-c1x;cut1y(1)-c1y];
vecarCrossWr = cross([vecar;0],[wr;0]);
if ( vecarCrossWr(3) > 0 )
   wr = -wr;
end

% Normal direction.
nwr = [wr(2) -wr(1)];
% unitNwr = nwr/norm(nwr);

% Force point.
fp = ([xCutCtr1;yCutCtr1]*area1+[xCutCtr2;yCutCtr2]*area2)/sumArea;

return