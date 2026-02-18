function [polygonX,polygonY,initialX,initialY,initialF,...
  polygonSide,polygonisfree,m,mi,wallNum]=DEM_initial_blender
% PURPOSE: 
% Set up the geometry of a blender:
%  1 circular outer wall
%  1 rotating propeller inside, a circle with 2 blades
%  Particles in the space between the propeller and the outer wall
%
% 1. Read in the geometry of the walls
% 2. Initialize the particles
% 3. Compute the masses, moments of inertia of the particles and wall
% 4. Concatenate the data for the respective particles and walls so that 
%    both can be treated together 
% OUTPUT: 
%  polygonX,polygonY contains the corner positions,
%  initialX,initialY,initialF contains the positions and their higher time derivatives
%  polygonSide is the number of corners per particle
%  polygonisfree is whether the particle is wall (0) or not(1)
%  m / mi are the mass and moment of inertia lists
% CAVEAT: 
%  - If two particles are initialized with overlaping bounding boxes, the 
%    neighborhood algorithm may fail to detect this pair of interacting particles
%  - Wall particles are considered to be "infinitely heavy" and are rescaled
%    to 1000 times the largest (physical) mass of any particle or wall
% PERFORMANCE CONSIDERATION: 
%  The MATLAB-editor gives warnings for several variables that the size
%  changes in every loop; because this happens only for the initialization,
%  no performance optimization was attempted
% TODO: 
%  Write an additional routine which makes sure that when a new particle is
%  initialized, it has no overlap with the previous ones; This can be done 
%  by the call of the overlap-routine for the newly initialized particle and
%  a loop over all previously initialized particles
% LITERATURE: sec. 7.4
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014
% ACTUAL DATA:
% r Drum:462.5mm
% r Wheel:307.5mm
% r Wing:437.5mm

global rhoPoly
global rhoWall
global gravity
global iblade1 iblade2
global r_i_i r_i_o
global nwall_i

% initialization of inner disk: 
% inner and outer radius, number of particles
BASE_R = 0.361;
nwall_i = 284; % MEMO: Depends on the size of drum 
iymax=36;

r_i_i = BASE_R-.007;
r_i_o = BASE_R;
rwheel = 0.4375;
iinner = 0;
for i = 1:nwall_i
  iinner = iinner + 1;
  numcorn = 4;%textscan(walldatain,'%d1','commentStyle','%');
  wallSide(i) = numcorn;
  wallisfree(i,1) = 0;
  iSector = i + 80; % FIXME: Added 80 to prevent rounding error
  alphastart = iSector * 2 * pi / nwall_i;
  alphaend = (iSector + 1) * 2 * pi / nwall_i;
  wallX(1,i)=r_i_i*cos(alphastart);
  wallY(1,i)=r_i_i*sin(alphastart);
  if (iinner==nwall_i/4-3)
    iblade1=i;
    wallX(2,i)=r_i_o*cos(alphastart);
    wallY(2,i)=r_i_o*sin(alphastart);
    wallX(3,i)=r_i_o*cos(alphaend);
    wallY(3,i)=r_i_o*sin(alphaend);
  elseif (iinner==3*nwall_i/4-3)
    iblade2=i;  
    wallX(2,i)=r_i_o*cos(alphastart);
    wallY(2,i)=r_i_o*sin(alphastart);
    wallX(3,i)=r_i_o*cos(alphaend);
    wallY(3,i)=r_i_o*sin(alphaend);
  else
    wallX(2,i)=r_i_o*cos(alphastart);
    wallY(2,i)=r_i_o*sin(alphastart);
    wallX(3,i)=r_i_o*cos(alphaend);
    wallY(3,i)=r_i_o*sin(alphaend);
  end
  wallX(4,i)=r_i_i*cos(alphaend);
  wallY(4,i)=r_i_i*sin(alphaend);
end  
wallNum = nwall_i;
disp(['Number of walls: ' num2str(wallNum)])

% xv, yv used inpolygon()
circleX = [wallX(1, 1:nwall_i) wallX(1, 1)];
circleY = [wallY(1, 1:nwall_i) wallY(1, 1)];

% a, b is diameter, not radius
a_min = 0.008;
a_max = 0.008;
b_min = 0.008;
b_max = 0.011;

max_d = max([a_min a_max b_min b_max]);
BUFFER = 1.1;

PARTICLE_CREATION_RATE = 258; % MEMO: Depends on the size of drum
xPosition = -PARTICLE_CREATION_RATE * max_d / 2 : max_d * BUFFER : PARTICLE_CREATION_RATE * max_d / 2;
yPosition = - (r_i_i - max_d / 2) : max_d * 0.9 : (r_i_i - max_d / 2);

% Initialize the maximal particle radius with 
% the width of half the left wall

particleNum = 0;
IX_MAX = 233; % MEMO: Depends on the size of drum
for iy = 1:iymax
  if 2*round(iy/2)==iy
    ixmax = IX_MAX + 1;
    dxoff = 0;
  else
    ixmax = IX_MAX;
    dxoff = 0.5 * max_d;
  end 
  for ix=1:1:ixmax   %space out particle centers in x
    particleNum=particleNum+1;
    ncorner = randi([7 12], 1); % randomize the number of corners
    % angles for the corners
    angles = 2 * pi / ncorner * (0:ncorner - 1);
    % Elongated particles, ellipse radius:
    ell_ar=(a_min+(a_max-a_min)*rand(1))/2;
    ell_br=(b_min+(b_max-b_min)*rand(1))/2;
    % rotate the particle around a random angle theta
    rand_ang = randi([0 179]);
    particleX(:,particleNum)=nan(12,1);
    particleY(:,particleNum)=nan(12,1);
    %  Center of mass
    for k = 1:ncorner
      theta = pi * 1 / 180 .* rand_ang;
      J = [cos(theta) -sin(theta); sin(theta) cos(theta)];
      XRel(k) = cos(angles(k))*ell_ar;
      YRel(k) = sin(angles(k))*ell_br;
      A = J * [XRel(k); YRel(k)];
      particleX(k,particleNum) = xPosition(ix) + A(1) + dxoff;
      particleY(k,particleNum) = yPosition(iy) + A(2);
    end
    particleSide(particleNum) = ncorner;
    particleisfree(particleNum, 1) = 1;
    % Initialize velocities as 0
    VX(particleNum, 1) = 0;
    VY(particleNum, 1) = 0;
    VF(particleNum, 1) = 0;

    In = inpolygon(particleX(1:ncorner,particleNum),...
      particleY(1:ncorner,particleNum), circleX, circleY);

    if sum(In)<ncorner
      particleX(:,particleNum)=0;
      particleY(:,particleNum)=0; 
      particleNum=particleNum-1 ;
    end
  end
end

if length(particleSide) > particleNum
  particleSide=particleSide(:,1:end-1);
  particleisfree=particleisfree(1:end-1,:);
end

disp(['Number of particles: ' num2str(particleNum)])

% Initialize masses, mom. of inertia, centers
% For Walls
[mwall,miwall,wallCOMF,wallCOMX,wallCOMY]=...
 mass_momentinertia(wallSide,wallX,wallY,rhoWall);

% For Particles
[m,mi,particleCOMF,particleCOMX,particleCOMY]=...
 mass_momentinertia(particleSide,...
 particleX,particleY,rhoPoly);
% Because the walls should be "infinitely heavy",
% their masses are set to "large" values (1000 times larger than the heaviest particle/wall)
maxm=max([m mwall]);
maxmi=max([mi miwall]);
for i=1:wallNum 
  mwall(i) =maxm;
  miwall(i)=maxmi;
end


% Combine data structures for particles and walls
maxSide = max([particleSide wallSide]); 
polygonNum = particleNum+wallNum;

disp(['Number of polygons (wall + particle): ' num2str(polygonNum)])

m=[m mwall]';
mi=[mi miwall]';
polygonX=zeros(maxSide,polygonNum);
polygonY=zeros(maxSide,polygonNum);
polygonSide=[wallSide particleSide];
polygonisfree=[wallisfree; particleisfree];

% Copy walls into the particle structure
for i=1:wallNum
  tempSide = wallSide(i);
  polygonX(1:tempSide,i) = wallX(1:tempSide,i);
  polygonY(1:tempSide,i) = wallY(1:tempSide,i);
end
iparticle=0;
% Copy particles into the particle structure
for i=wallNum+1:wallNum+particleNum
  iparticle=iparticle+1;
  tempSide=particleSide(iparticle);
  polygonX(1:tempSide,i)=particleX(1:tempSide,iparticle);
  polygonY(1:tempSide,i)=particleY(1:tempSide,iparticle);
end


% Variables for predictor the predictor corrector:
% up to six values, the positions and its 5 time derivatives:
% 1. Center of mass
% 2. Velocity
% 3. Acceleration
% 4. Jerk
% 5. First derivitive of jerk.
% 6. Second derivitive of jerk.
initialX = zeros(particleNum,6);
initialY = zeros(particleNum,6);
initialF = zeros(particleNum,6);
for i=1:wallNum
% wall velocities and their time derivatives are 0, if walls moving at 
% constant velocity are desired, changes must also be made in the predictor
% and the corrector, where in the current version the wall velocities are 
% set to zero.
  initialX(i,:)=[wallCOMX(i) 0 0 0 0 0]; % x-position
  initialY(i,:)=[wallCOMY(i) 0 0 0 0 0]; % y-position
  initialF(i,:)=[wallCOMF(i) 0 0 0 0 0]; % orientation
end
for i=1:particleNum
  initialX(wallNum+i,:)=[particleCOMX(i) VX(i) 0 0 0 0];        
  initialY(wallNum+i,:)=[particleCOMY(i) VY(i) -gravity 0 0 0]; 
  initialF(wallNum+i,:)=[particleCOMF(i) VF(i) 0 0 0 0];        
end



return