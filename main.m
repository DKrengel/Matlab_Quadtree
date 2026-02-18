% PURPOSE: Determine effect of blade shape for Blenders 
% NOTE:  Neighborhood is computed with a quadtree
% CAVEATS: No inclusion of drag, particle size not realistic
% AUTHORS: Shi Han Ng, Hans-Georg Matuttis,
%          Ko Kandori, Yuki Watanabe
% REFERENCES: HG Matuttis, J. Chen "Understanding the Discrete Element Method", Wiley
%             DOI:10.1002/9781118567210
% REVISION HISTORY:

clear all
clear neighborhoodAlgorithmTreeCode % Clear persistent variables in neighborhood_algorithm.
clear neighborhood_algorithm % Clear persistent variables in neighborhood_algorithm.
format compact

startingTime = tic; % cputimestart = cputime;

global gravity rhoPoly rhoWall emodul mu gamma  
global nwall_i omega_blender
global t
global iblade1 iblade2
global it


omega_blender=2*pi*0.25 % f= [1/s], rotation speed of the blender, converted into radian

% set seed to obtain the same configuration 
% in every program run: only for debugging
rand('seed',1) 

% Time-variables
t     = 0.0;            % Initial simulation time
tmax  = 0.20;           % Final simulation time 
dt    = 2.5e-5;          % Size of the timestep 
nt    =  round((tmax-t)/dt); % Number of steps
everyn=2500;           % Number of steps between graphical output

% DEM-parameters
gravity   = 9.81;    % [m/s^2] Standard acceleration due to free fall.
rhoPoly   = 2500.0;  % [kg/m^3] Density of particle.
rhoWall   = rhoPoly; % [kg/m^3] Density of wall.
emodul    = 1.0e+7;  % Young modulus. javascript:dPlayNext();
mu        = 0.23;     % Coefficient of friction.
gamma     = .0251;     % Damping constant.

%--------------------------------------------------------------------------
% Initialization of the particles; initial coordinates named "corrected" to
% avoid recopying with other variable names for the main loop
[particleX,particleY,correctedX,correctedY,correctedF,...
 particleSide,freeParticle,m,mi,numWall] = DEM_initial_blender;

% Maximum stepsize according to the natural frequency of the lightest particle, p.240.
% Analogous to the discussion on 
% p. 291 for the three-dimensional case
MAXDT = 0.1*sqrt( min(m)/emodul )*pi
if ( dt > MAXDT )
  error('The stepsize is too large.'); 
elseif(dt<5*MAXDT)
  disp(['The timestep of ' num2str(dt) ])  
  disp('is considerably smaller than the')
  disp(['maximally allowed timestep of ' num2str(MAXDT) ])  
end

% Graphical output of the initial configuration,
% to make sure that particles are where they should be
figure(1), clf, hold on
draw_particle(particleSide,particleX,particleY,numWall);
axis equal

disp('press key to continue')
pause


%--------------------------------------------------------------------------


% Main loop
disp('Starting time integration')
for it=1:nt

% Increment the simulation time.
  t=t+dt;

% Predictor step 
  [PX,VX,PY,VY,PF,VF] = gearPredictor(dt,freeParticle,...
   correctedX,correctedY,correctedF);

% Update the vertices according to the predictor coordinates
  [particleX,particleY,PX,PY,VX,VY,VF] = updatevertices_blender(...
    particleSide,correctedX,correctedY,correctedF,...
    PX,PY,PF,VX,VY,VF,particleX,particleY,dt);

 % Neighborhood algorithm based on sort-and-sweep
 % [col_list, col_list_len] = neighborhood_algorithm(...
 %      particleSide, particleX, particleY);

 % Neighborhood algorithm based on the tree-code
  [col_list, col_list_len] = neighborhoodAlgorithm_QTree(...
       particleSide, particleX, particleY, numWall);

 % Interaction computation 
  [acx,acy,acf,overlapList]=DEM_interaction(col_list,col_list_len,...
        dt,particleX,particleY,particleSide,...
        PX,PY,VX,VY,VF,freeParticle,m,mi);

% Corrector-step of the Gear Predictor-Corrector 
  [correctedX,correctedY,correctedF]=...
  gearCorrector(acx,acy,acf,dt,freeParticle);

% Graphics output every everyn steps
  if (everyn*round(it/everyn)==it)    
    it
  %  figure(2), clf, hold on
  %    draw_particle(particleSide,particleX,particleY,numWall)
  %   axis image
  %    drawnow
  end

end % end of time evolution loop


disp('Finished without error!')
%Display the final configuration
figure(2), clf, hold on
draw_particle(particleSide,particleX,particleY,numWall);
axis equal

return