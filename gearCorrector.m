function [XC,YC,FC]=gearCorrector(AX,AY,AF,dt,freeParticle)
% PURPOSE:
%  corrector-part of the Gear Predictor-Corrector,
%        | xc|   | xp|     |c0|
%        |vxc|   |vxp|     |c1|
%        |axc| = |axp|  +  |c2|.
%        |bxc|   |bxp|     |c3|
%        |cxc|   |cxp|     |c4|
%        |dxc|   |dxp|     |c5|
% USAGE: after the force calculation;
%  comment out either the fith order or the second order coefficients to change the order
% INPUT:
%  Corrected values of the second time derivatives of the coordinates  
%  (acceleration/forces rescaled with masses).
%  AX,AY,AF (x-coordinate,y-coordinate, orientation)
%  dt - Timestep.
%  freeParticle is 0 or 1 depending of whether the particle moves freely, 
%  or not (in which case  there are no position changes)
%  XP,YP,FP are the predicted positions and their time derivatives, passed 
%  as global variables
% OUTPUTS: The "un-scaled" (SI-units) corrected values 
%  XC,YC,FC of the position/orientation and their time derivatives.
% CAVEAT: Higher order does not necessarily mean higher accuracy, but higher stability. 
% RECORD OF REVISION:
%  13-May-2014   HG Matuttis   inclusion of isfree
%  2013-May-30   S. H. Ng      Starts to record.
% LITERATURE: See comments for gearPredictor

global nwall_i
global XP YP FP

% set the accerlations to zero for wall particles
AX = AX.*freeParticle;
AY = AY.*freeParticle;
AF = AF.*freeParticle;

% Uncomment this and comment out the next assignment of C if you want to 
% work with the 5th order method;  don't forget to change the corrector in 
% gearPredictor.m, too.
% Coefficients for 5th-order Predictor-Corrector
 C0 = 3.0   / 16.0;
 C1 = 251.0 / 360.0;
 C2 = 1.0;
 C3 = 11.0  / 18.0;
 C4 = 1.0  / 6.0;
 C5 = 1.0  / 60.0;

% Coefficients for 2nd-order Predictor-Corrector
%C0 = 0.0;
%C1 = 1.0;
%C2 = 1.0;
%C3 = 0.0;
%C4 = 0.0;
%C5 = 0.0;

pdt2 = 0.5 * dt * dt;

idt = 1/dt;
ipdt2 = idt * 2.0 / dt;
ipdt3 = ipdt2 * 3.0 / dt;
ipdt4 = ipdt3 * 4.0 / dt;
ipdt5 = ipdt4 * 5.0 / dt;

corrCoeff = [C0;C1;C2;C3;C4;C5];

% Scale AX, AY and AF and compute the deltas.
deltaX = AX'*pdt2 - XP(3,:);
deltaY = AY'*pdt2 - YP(3,:);
deltaF = AF'*pdt2 - FP(3,:);

%size(XC)
%size(XP)
%size(DeltaX)
XC=zeros(size(XP));
YC=zeros(size(YP));
FC=zeros(size(FP));


XC(:,nwall_i+1:end) = XP(:,nwall_i+1:end) + corrCoeff*deltaX(nwall_i+1:end);
YC(:,nwall_i+1:end) = YP(:,nwall_i+1:end) + corrCoeff*deltaY(nwall_i+1:end);
FC(:,nwall_i+1:end) = FP(:,nwall_i+1:end) + corrCoeff*deltaF(nwall_i+1:end);
XC(:,1:nwall_i) = XP(:,1:nwall_i);
YC(:,1:nwall_i) = YP(:,1:nwall_i);
FC(:,1:nwall_i) = FP(:,1:nwall_i);


% "Un-scaled" corrected values.
XC(2,nwall_i+1:end) = XC(2,nwall_i+1:end)*idt;
XC(3,nwall_i+1:end) = XC(3,nwall_i+1:end)*ipdt2;
XC(4,nwall_i+1:end) = XC(4,nwall_i+1:end)*ipdt3;
XC(5,nwall_i+1:end) = XC(5,nwall_i+1:end)*ipdt4;
XC(6,nwall_i+1:end) = XC(6,nwall_i+1:end)*ipdt5;

YC(2,nwall_i+1:end) = YC(2,nwall_i+1:end)*idt;
YC(3,nwall_i+1:end) = YC(3,nwall_i+1:end)*ipdt2;
YC(4,nwall_i+1:end) = YC(4,nwall_i+1:end)*ipdt3;
YC(5,nwall_i+1:end) = YC(5,nwall_i+1:end)*ipdt4;
YC(6,nwall_i+1:end) = YC(6,nwall_i+1:end)*ipdt5;

FC(2,nwall_i+1:end) = FC(2,nwall_i+1:end)*idt;
FC(3,nwall_i+1:end) = FC(3,nwall_i+1:end)*ipdt2;
FC(4,nwall_i+1:end) = FC(4,nwall_i+1:end)*ipdt3;
FC(5,nwall_i+1:end) = FC(5,nwall_i+1:end)*ipdt4;
FC(6,nwall_i+1:end) = FC(6,nwall_i+1:end)*ipdt5;

XC = XC';
YC = YC';
FC = FC';

return

