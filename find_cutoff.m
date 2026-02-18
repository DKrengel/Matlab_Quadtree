function [sxcut, sycut, sarea, subpolygonx, subpolygony, xCutoffCenter, ...
    yCutoffCenter] = find_cutoff(sx, sy, X1, Y1, Xc, Yc)
%  find_cutoff
%  PURPOSE:
%     Output those two points  which are "cut off"
%     by the other polygon; find_cutoff has to
%     be called independent for the two contacting polygons
%
%     Input the information of a polygon
%     with the vertex coordinates X1,Y1
%     and the coordinates of the centroid X1c,Y1c
%     as well as several points sx(i),sy(i)
%     sxcut,sycut will those points from sx(i),sy(i)
%     between which all other points of sx(i),sy(i)
%     are situated as seen from the centroid X1c,Y1c
%     sarea is the area between the points sx(i),sy(i)
%
%  USAGE:
%     Call if the overlap-routine gives more than two intersection 
%     point between polygons
%
%  RECORD OF REVISION:
%         Date       Programmer      Description of change
%         ====       ==========      =====================
%      2013-Apr-22   H.-G. Matuttis  Original.
%   1. 2013-Apr-22   S. H. Ng        Bug fix: preserve the last entry of 
%                                    "corneralpha" because there is no
%                                    periodicity in the elements of
%                                    "corneralpha".
%   2. 2016-07-27    Masaki Nawa     Bug fix: transpose the array of 
%                                    "subpolygonx" and "subpolygony".[Rev:2]
%
%  ALGORITHM:
%     1. Compute the auxiliary point sx_a,sy_a as
%        the average of the sx,sy
%     2. Compute the rays from X1c,X2c to the (sx(i),sy(i))
%     3. Compute the angles of these rays with
%        respect to (sxc,syc) to the (sx(i),sy(i))
%        via atan2(py,px), i.e. via 
%        the inner product between (sxc,syc) to the (sx(i),sy(i)) (px) 
%        and the inner product between the normal of (sxc,syc)=(-syc sxc)
%        and  the (sx(i),sy(i)) (py) 
%        The computation of the atan2 with respect to 
%        sx_a-Xc, sy_a-Yc is necessary to make sure that the 
%        angles are in the right quadrants
%     4. Choose those sx(i),sy(i) with of the rays with the extremal 
%        angles
%
%  CAVEATS: 
%     If the centroid is inside the polygon defined by sx,sy, 
%     the algorithm will fail!
%

coder.varsize('subpolygonx');
coder.varsize('subpolygony');

if (length(sx)~=length(sy))
  error('sx and sy must have the same length')  
end 

if (length(X1)~=length(Y1))
  error('X1 and Y1 must have the same length')  
end    

if (length(sx)<=2)
  error('sx must have at least 3 elements ')
end

% Comment out "if ( shouldOutputFig )"-block for production run.
% shouldOutputFig = 1;
% if ( shouldOutputFig )
%    figure;
%    clf
%    hold on
%    hParticleOutline = plot([X1;X1(1)], [Y1;Y1(1)], '-k');
%    set(get(get(hParticleOutline,'Annotation'), 'LegendInformation'), ...
%     'IconDisplayStyle', 'off') % Exclude outline from legend.
%    plot(Xc, Yc, 'xk', 'MarkerSize', 12, 'DisplayName', 'Xc/Yc')
%    plot(sx, sy, 'xb', 'MarkerSize', 12, 'DisplayName', 'Intersection points');
%    hold off
% end

sx_a=mean(sx);
sy_a=mean(sy);

% if ( shouldOutputFig )
%    hold on
%    plot(sx_a, sy_a, '*', 'DisplayName', 'sx_a/sy_a')
%    hold off
% end

% Ray from the centroid of the polygon to the average
% of the intersection points
% cax=(sx_a-Yc);
% cay=(sy_a-Xc);
cax=(sx_a-Xc);
cay=(sy_a-Yc);
% if ( shouldOutputFig )
%    hold on
%    quiver(Xc, Yc, cax, cay, 'AutoScale', 'off', 'DisplayName', 'cax/y');
%    hold off
% end

% % Test: For the line between the centroid of the polygon
% %       and the average of the intersection points, the angle
% %       should be zero
% plot([sx_a Xc],[sy_a Yc],'k--')
% px=dot([cax cay],[cax cay])                                   
% py=dot([cax cay],[-cay cax])                                 
% alpha=atan2(py,px)
% text(sx_a,sy_a,num2str(alpha/pi*180),'Fontsize',18)

% Compute the angles for the rays
sxLen = length(sx);
rx    = zeros(1, sxLen);
ry    = zeros(1, sxLen);
alpha = zeros(1, sxLen);
for i=1:sxLen
   rx(i)=(sx(i)-Xc);
   ry(i)=(sy(i)-Yc);
   px=dot([rx(i) ry(i)],[cax cay]);
   py=dot([rx(i) ry(i)],[-cay cax]);
   alpha(i)=atan2(py,px);
   % % Plot the respective rays, with the respective angles
   % % (in degrees) to the line between the centroid of the polygon
   % % and the average of the intersection points
   % hold on
   % plot([sx(i) Xc],[sy(i) Yc],'k-')
   % text(sx(i),sy(i),num2str(alpha(i)/pi*180),'Fontsize',18)
   % hold off
end
% if ( shouldOutputFig )
%    hold on
%    quiver(repmat(Xc,size(rx)), repmat(Yc,size(rx)), rx, ry, ...
%       'AutoScale', 'off', 'DisplayName', 'rx/y');
%    text(sx, sy, num2str(alpha'/pi*180), 'VerticalAlignment', 'top', ...
%       'FontSize', 12)
%    hold off
% end


% Look for the rays with the largest and smallest angle.
% #################################
% maxangle=alpha(1);
% maxangle_index=1;
% minangle=alpha(1);
% minangle_index=1;
% for i=2:length(sx)
%   if (alpha(i)>maxangle)
%     maxangle=alpha(i);
%     maxangle_index=i;
%   end
%   if (alpha(i)<minangle)
%     minangle=alpha(i);
%     minangle_index=i;
%   end                                                                                      
% end  
% ################################# 
% The above is replaced with the following:
[~, maxangle_index] = max(alpha);
[~, minangle_index] = min(alpha);
% plot(sx(maxangle_index),sy(maxangle_index),'p','Markersize',16)
% plot(sx(minangle_index),sy(minangle_index),'p','Markersize',16)
sxcut=zeros(1,2);
sycut=zeros(1,2);
sxcut(1)=sx(minangle_index);
sycut(1)=sy(minangle_index);
sxcut(2)=sx(maxangle_index);
sycut(2)=sy(maxangle_index);

% if ( shouldOutputFig )
%    hold on
%    plot(sxcut, sycut, 'og', 'MarkerSize', 15, 'DisplayName', 'Cutoff points')
%    hold off
% end

% % Compute the area of the cut off polygon:
% % Sum up the area of all the triangles
% % from the averages of the (sx(i),sy(i))
% % to the (sx(i),sy(i))
% spx=sx;
% spx(end+1)=sx(1);
% spy=sy;
% spy(end+1)=sy(1);
% text(spx(1,:),spy(1,:),num2str([1:length(spx(1,:))]'))
% text(spx(2,:),spy(2,:),num2str([1:length(spx(2,:))]'))

% Compute  total area of the particle
particle_area=abs(oriented_polygon_area(X1,Y1));

% Find out which corners of the polygon are between sxcut and sycut
corneralpha=zeros(size(X1));
for i=1:length(X1)
   xray=X1(i)-Xc;
   yray=Y1(i)-Yc;
   corneralpha(i)=atan2(yray,xray);
   if (corneralpha(i)<0)
      corneralpha(i)=corneralpha(i)+2*pi;
   end
   % % Plot the respective rays, with the respective angles
   % % (in degrees) to the line between the centroid of the polygon
   % % and the average of the intersection points
   % plot([X1(i) Xc],[Y1(i) Yc],'k-')
   % text(X1(i),Y1(i),[num2str(i) '/' num2str(corneralpha(i)/pi*180,3)],'Fontsize',18)
end
% if ( shouldOutputFig )
%    hold on
%    corneralphaStr = arrayfun(@(i) ...
%       {sprintf('%i/%.2f', i, corneralpha(i)/pi*180)}, (1:length(X1))');
%    text(X1, Y1, corneralphaStr,'Fontsize',18)
%    hold off
% end

intersalpha=zeros(1,2);

xray=sxcut(1)-Xc;
yray=sycut(1)-Yc;
intersalpha(1)=atan2(yray,xray);
if (intersalpha(1)<0)
  intersalpha(1)=intersalpha(1)+2*pi;
end
xray=sxcut(2)-Xc;
yray=sycut(2)-Yc;
intersalpha(2)=atan2(yray,xray);
if (intersalpha(2)<0)
  intersalpha(2)=intersalpha(2)+2*pi;
end

% if ( shouldOutputFig )
%    hold on
%    hl1 = plot([sxcut(1) Xc],[sycut(1) Yc],'k-');
%    hl2 = plot([sxcut(2) Xc],[sycut(2) Yc],'k-');
%    set(get(get(hl1,'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
%    set(get(get(hl2,'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
%    text(0.5*(sxcut(1)+Xc),0.5*(sycut(1)+Yc),sprintf('sc1-%.2f', intersalpha(1)/pi*180),'Fontsize',18)
%    text(0.5*(sxcut(2)+Xc),0.5*(sycut(2)+Yc),sprintf('sc2-%.2f', intersalpha(2)/pi*180),'Fontsize',18)
%    hold off
% end

% Find the groups of corneralpha-values which are 
% limited by the intersalpha-values
% Bring the angles in an order from 0 to 2pi
% corneralpha=corneralpha(1:end-1); % Commented out because there is no
                                    % repetition in the last entry.
[cornersort,cornerind]=sort(corneralpha);
cornersort(end+1)=cornersort(1)+2*pi;

% Exchange the intersection points if necessary
% so that they are ordered respective to the angle
if ( intersalpha(1) > intersalpha(2) )
   intersdum=intersalpha(2);
   intersalpha(2)=intersalpha(1);
   intersalpha(1)=intersdum;
   sxdum=sxcut(2);
   sxcut(2)=sxcut(1);
   sxcut(1)=sxdum;
   sydum=sycut(2);
   sycut(2)=sycut(1);
   sycut(1)=sydum;
end

startgroup1=-1;
if ( intersalpha(1) < cornersort(1) )
   startgroup1=0;
end

if ( startgroup1 > -1 )
   startgroupset=1;
else
   startgroupset=0;
end

endgroup1=0;%dummy:definitely difined later 
for i=1:length(cornersort)-1
   % [cornersort(i) cornersort(i+1)];
   if startgroupset==0
      if ( intersalpha(1)>cornersort(i))&&(intersalpha(1)<cornersort(i+1) )
         startgroup1=i;
         startgroupset=1;
      end
   else
      if ( intersalpha(2)>cornersort(i))&&(intersalpha(2)<cornersort(i+1) )
         endgroup1=i;
         break
      end
   end
end
group1index=startgroup1+1:endgroup1;
% cornersort(group1index)*180/pi;

subpolygonx=sxcut(1);
subpolygony=sycut(1);
n=1;
for i=1:length(group1index)
   n=n+1;
   subpolygonx(n)=X1(cornerind(group1index(i)));
   subpolygony(n)=Y1(cornerind(group1index(i)));
end
n=n+1;
subpolygonx(n)=sxcut(2);
subpolygony(n)=sycut(2);
cutoffarea=abs(oriented_polygon_area(subpolygonx,subpolygony));
if ( cutoffarea > particle_area )
   disp('##### Ops! #####')
   %disp('Cutoff area cannot be larger than particle area.')
   %disp('Entering debug-mode >>>')
   error('Cutoff area cannot be larger than particle area.')%exit
end   

% if ( shouldOutputFig )
%    hold on
%    patch(subpolygonx, subpolygony, 1, 'FaceAlpha', 0.5, ...
%       'DisplayName', 'Cutoff polygon')
%    hold off
% end

% choose that area which is cut of by the sx(i),sy(i) wich is smaller and which presumably does not contain the center of mass
[xSubPolygonCenter, ySubPolygonCenter, areaSubPolygon] = centerOfMass(...
   subpolygonx, subpolygony);
if ( abs(cutoffarea-areaSubPolygon) > 1.0e-16 )
   disp('##### Ops! #####')
   error('Computation of the area of subpolygon is not consistance.')
   %disp('Entering debug-mode >>>')
   %exit
end

difference_area=particle_area-cutoffarea;
if (difference_area<cutoffarea)
   % Current subpolygon is inverted part of the cutoff.
   sarea=difference_area;
   % xCutoffCenter = (Xc*particle_area-xSubPolygonCenter*cutoffarea)/difference_area;
   % yCutoffCenter = (Yc*particle_area-ySubPolygonCenter*cutoffarea)/difference_area;
   
   % Attemp to compute the correct geometry of the cutoff.
   periodicCornerInd = [cornerind; cornerind];
   group1InvertedId = (endgroup1+1):(startgroup1+size(cornerind,1));
   cutCornerX1 = X1(periodicCornerInd(group1InvertedId));
   cutCornerY1 = Y1(periodicCornerInd(group1InvertedId));
   subpolygonx = [sxcut(2); cutCornerX1; sxcut(1)];
   subpolygony = [sycut(2); cutCornerY1; sycut(1)];
%    if ( shouldOutputFig )
%       hold on
%       patch(subpolygonx, subpolygony, 1, 'FaceAlpha', 0.5, ...
%          'FaceColor', 'y', 'DisplayName', 'Cutoff polygon')
%       hold off
%    end
   [xSubPolygonCenter, ySubPolygonCenter, areaSubPolygon] = centerOfMass(...
      subpolygonx', subpolygony');
   if ( abs(sarea-areaSubPolygon) > 1.0e-16 )
      disp('##### Ops! #####')
      error('The computation of the cutoff polygon might be wrong.')
      %disp('Entering debug-mode >>>')
      %exit
   end
   xCutoffCenter = xSubPolygonCenter;
   yCutoffCenter = ySubPolygonCenter;
   
else
   sarea=cutoffarea;
   xCutoffCenter = xSubPolygonCenter;
   yCutoffCenter = ySubPolygonCenter;
end

% if ( shouldOutputFig )
%    hLegend = legend('show');
%    set(hLegend, 'Interpreter', 'none');
%    axis square
% end

return

