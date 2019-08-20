%% Init

clear
close all
clc

addpath(genpath('./'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change flight parameters

coaxial.state.axial_vel = 0; %m/s 
coaxial.state.forward_vel = 0; %m/s 
coaxial.state.side_vel = 0; %m/s

% Control Inputs

coaxial.state.collective_u = 20; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.collective_l = 20; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.cyclic_s = 0; %sine term for cyclic (gets multiplied by sin(azimuth))
coaxial.state.cyclic_c = 0; %cosine term for cyclic (gets multiplied by cos(azimuth))

epsilon = 0.0001; %convergence accuracy for Prandtl tip function and inflow ratio

%warning('off')

% Testing 

plots= true;
verbose= true;
debug = false;
method='leishman'; %'leishman','airfoil'

[Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);

%legend below (SWIER)
%{
%Power is a 1x2 matrix -> [P_upper; P_lower] [W]

%Forces is a 3x2 matrix -> [Fx_upper, Fy_upper, Fz_upper; Fx_lower, Fy_lower, Fz_lower] [N]
%using the right-handed coordinate system with x pointing forward, z pointing up

%Moments is a 3x2 matrix -> [Mx_upper, My_upper, Mz_upper; Mx_lower, My_lower, Mz_lower] [Nm] 
%using the same coordinate system as Forces

%CT is the thrust coefficient 1x2 matrix -> [CT_u; CT_l];
%CP is the torque coefficient 1x2 matrix -> [CP_u;CP_l];

%} 
