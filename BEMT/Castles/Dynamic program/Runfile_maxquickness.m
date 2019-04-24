% Runfile for dynamic pitch response of the coaxial helicopter
clear all
clc
close all

% Loading trim conditions
load coax_trimstates
Ka32data
m2kts=1852/3600;
% Defining intitial speed
Vin=20;

% Combining intial speed to tabulated data
N=round(Vin)+1;

in_states=states(N,:)';
in_controls=controls(N,:)';

xInitial=[in_states];

% Defining required attitudes
thetaFreq=30*pi/180+in_states(8);
t_theta=0.2;

% ureq=70*m2kts;
% t_u=0.2;
% t_u1=0;
% ureq1=0*m2kts;
% t_u2=0;
% ureq2=70*m2kts;
% t_u3=30;
% ureq3=20*m2kts;
% t_u4=60;
% ureq4=20*m2kts;


phiFreq=in_states(9);
t_phi=0;

hreq=-in_states(12);
zreq=-hreq;
t_z=0;

%t_rreq=0.2;
%rreq=0.01;
hdreq=in_states(7);
t_hd=0;

% Defining gains

% Longitudinal cyclic
Ktheta=-3.1;
Ktheta_int=-0.1;
Kq=0.6;
% Ku=-0.1;
% Ku_int=-0.0001;


% Lateral cyclic
Kphi=0.6;
Kphi_int=-0.02;
Kp=-0.2;

% Collective thrust

Kz=-0.005;
Kz_int=0;
Kzdot=0;%-0.01;

% Yaw control
Kpsi=0.7;
Kpsi_int=-0.001;
Kr=-0.4;%0.01;
Kdiffu=0.5;
Kdiffl=0.5;


[t,x,y]=sim('PIDcoaxial',[0 6])





