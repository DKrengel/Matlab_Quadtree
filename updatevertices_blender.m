function [particleX,particleY,PX,PY,VX,VY,VF]=...
  updatevertices_blender(particleSide,...
  PXOld,PYOld,PFOld,PX,PY,PF,VX,VY,VF,particleX,particleY,dt);
% PURPOSE: 
%  Rotation with respect to the center of mass and translation according to 
%  the difference between the current and the previous timestep
% USAGE: 
%  If a predictor-corrector integrator is used, in timestep n it must be 
%  called after the predictor-step, as the corrected coordinates of 
%  timestep (n-1) and the predicted coordinates of timestep n are used.
% ALGORITHM:
%  1. Shift the center of mass into the origin,
%  2. Rotate the particle with an angle which is the difference between the
%     old and the new orientation
%  3. Translate the corners to the new center of mass (predicted coordinate)
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

global nwall_i omega_blender
global XP YP
global iblade1 iblade2
global it


particleNum=length(PX);
particleRelX=zeros(size(particleX));
particleRelY=zeros(size(particleY));

for i=nwall_i+1:particleNum
% Relative Position=Center of Mass - Absolute Position
  particleRelX(1:particleSide(i),i) = PXOld(i) - particleX(1:particleSide(i),i);
  particleRelY(1:particleSide(i),i) = PYOld(i) - particleY(1:particleSide(i),i);
  ps=particleSide(i);
  % ROTATION AND PARTICLES UPDATES
  [particleX(1:ps,i),particleY(1:ps,i)]=...
  rotation(particleRelX(1:ps,i),particleRelY(1:ps,i),...
  PX(i,1),PY(i,1),PF(i,1)-PFOld(i));
end

% The wall(drum) is not turned until particles are deposited.
if it < 30000
  return
end

Arotate=[cos(omega_blender*dt) -sin(omega_blender*dt)
         sin(omega_blender*dt)  cos(omega_blender*dt)];
for i = 1:nwall_i
  r1=[PXOld(i)
      PYOld(i)];    
  r2=Arotate*r1;
  PX(i,1)=r2(1);
  PY(i,1)=r2(2);
  XP(1,i)=r2(1);
  YP(1,i)=r2(2);

  particleXY=[ particleX(1:particleSide(i),i)'
               particleY(1:particleSide(i),i)'];

  rotatedvec=(Arotate*particleXY)';
  particleX(1:particleSide(i),i)=rotatedvec(:,1);
  particleY(1:particleSide(i),i)=rotatedvec(:,2);
 
end

VX(1:nwall_i)=-PY(1:nwall_i,1)*omega_blender;
VY(1:nwall_i)= PX(1:nwall_i,1)*omega_blender;
VF(1:nwall_i)=omega_blender;

return
