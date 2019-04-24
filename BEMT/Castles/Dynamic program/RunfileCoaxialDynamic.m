      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %    Non-linear, 6DOF dynamic simulation for the Coaxial helicopter %
      %                                                                   %
      %                     Master of science thesis                      %
      %                         February 2008                             %
      %                                                                   %
      %                     By John Campfens 1296434                      %
      %                 Delft University of Technology                    %
      %                 Faculty of Aerospace Engineering                  %
      %                                                                   %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Runfile for dynamic pitch response of the coaxial helicopter
clear all
clc
close all

% Loading trim conditions
load coax_trimstates
Ka32data

% Defining intitial speed
Vin=20;

% Combining intial speed to tabulated data
N=Vin+1;

in_states=states(N,:)';
in_controls=controls(N,:)';

xInitial=[in_states];

% Defining required attitudes
thetaFreq=in_states(8);
t_theta=0;

phiFreq=in_states(9);
t_phi=0;

hreq=-in_states(12);
zreq=-hreq;
t_z=0;

hdreq=in_states(7);
t_hd=0;

% Defining gains

% Longitudinal cyclic
Ktheta=-0.8;
Ktheta_int=-0.1;
Kq=0.3;

% Lateral cyclic
Kphi=0.4;
Kphi_int=0.002;
Kp=-0.3;

% Collective thrust

Kz=-0.005;
Kz_int=0;
Kzdot=0;

% Yaw control
Kpsi=0.7;
Kpsi_int=-0.001;
Kr=-0.4;
Kdiffu=0.5;
Kdiffl=0.5;


[t,x,y]=sim('PIDcoaxial',[0 5])





