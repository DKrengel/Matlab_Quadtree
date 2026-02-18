function draw_particle(s, x, y, wallCount)
%  PURPOSE: Plot the polygonal particles and walls
%  REVISION HISTORY: S. H. Ng, HG Matuttis 09-May-2014
%  TODO: Plot particles and walls in different color

global t
n=length(x);

hold off
hold on       
for i = 1:n
  ncorner = s(i);
  if i == 1
    patch(x(1:ncorner,i),y(1:ncorner,i),'b');
  else
    patch(x(1:ncorner,i),y(1:ncorner,i),'g');
  end
% remove the comment to annotate the corner number 
% for each particle
%  text(x(1:ncorner,i),y(1:ncorner,i),num2str([1:ncorner]'));
%  text(mean(x(1:ncorner,i)),mean(y(1:ncorner,i)),num2str(i));
end
% for i = 1:wallCount
%   ncorner = s(i);
%   text(mean(x(1:ncorner, i)), mean(y(1:ncorner, i)), num2str(i), 'Color', 'w', 'HorizontalAlignment', 'center')
% end
axis equal
%drawnow
%eval(['print -depsc2 blender' num2str(t) '.eps'])

return