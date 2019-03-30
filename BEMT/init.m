%This file should be run before anything else in the pipeline. It puts all
%the necessary parameters in the workspace.

clear;
clc;
close all;

%% Atmospheric parameters

atm.rho = 1.225; %kg/m^3 sea-level air density
atm.gamma_air = 1.4; %[-] for air (diatomic gas value)
atm.R_air = 8.314510/0.0289645; %J/kg/K= N m / N s^2 m^-1 /K = m^2/s^2/K
atm.Temp = 20; %?C 
atm.c0 = sqrt(atm.gamma_air*atm.R_air*(atm.Temp+273.15)); %m/s - speed of sound

%% Coaxial rotor parameters

% Top rotor
rotor(1).cl_alpha = 2*pi;         %1/rad - Lift slope, thin airfoil theory for now
rotor(1).Nb = 2;             %# Number of blades
rotor(1).diameter = 3.6;             %m - using reference values from papers for now
rotor(1).rpm = 800;             %rev/min - using reference values from papers for now
rotor(1).chord = 0.15;             %m - using reference values from papers for now
rotor(1).solidity = rotor(1).Nb*rotor(1).chord/(pi*rotor(1).diameter*0.5);   % Blade solidity - not sure what to take as chord for tapered blade

% Bottom rotor
rotor(2).cl_alpha = rotor(1).cl_alpha;  %1/rad - Lift slope, thin airfoil theory for now
rotor(2).Nb = 2;   %# Number of blades
rotor(2).diameter = 3.6;             %m - using reference values from papers for now
rotor(2).rpm = 800;             %rev/min - using reference values from papers for now
rotor(2).chord = 0.15;             %m - using reference values from papers for now
rotor(2).solidity = rotor(2).Nb*rotor(2).chord/(pi*rotor(2).diameter*0.5);   % Blade solidity - not sure what to take as chord for tapered blade
rotor(2).rd = 0.82*rotor(2).diameter/2;  % annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
rotor(2).A_Ad = 1/rotor(2).rd^2;


% Drag coefficients
Cd0 = 0.011;
D1 = 0.01;
D2 = 0.01;


%% Testing

%Prandtl_tip_loss(0.5,0,rotor(1))
