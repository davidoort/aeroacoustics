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

%% Top rotor
rotor(1).Nb = 2;             %# Number of aeros
rotor(1).R = 1.8;             %m - using reference values from papers for now
rotor(1).rpm = 800;             %rev/min - using reference values from papers for now
rotor(1).chord = 0.15;             %m - using reference values from papers for now
rotor(1).solidity = rotor(1).Nb*rotor(1).chord/(pi*rotor(1).R);   % aero solidity - not sure what to take as chord for tapered aero
% Aerodynamics
aero.cl_alpha = 2*pi;         %1/rad - Lift slope, thin airfoil theory for now
aero.Cd0 =  0.011;
aero.D1 = 0.01;
aero.D2 = 0.01;
rotor(1).aero = aero;

%Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
%CL_alpha
rotor(1).solidity = 0.027; 
aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
aero.Cd0 =  0.011;
aero.D1 = 0.0;
aero.D2 = 0.0;
rotor(1).aero = aero;


%% Bottom rotor
rotor(2).Nb = 2;  %# Number of aeros
rotor(2).R = 1.8;  %m - using reference values from papers for now
rotor(2).rpm = 800;  %rev/min - using reference values from papers for now
rotor(2).chord = 0.15;  %m - using reference values from papers for now
rotor(2).solidity = rotor(2).Nb*rotor(2).chord/(pi*rotor(2).R);   % aero solidity - not sure what to take as chord for tapered aero
rotor(2).rd = 0.82;  %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"

% Aerodynamics - maybe a bit of overkill at the moment
aero.cl_alpha = 2*pi;         %1/rad - Lift slope, thin airfoil theory for now
aero.Cd0 =  0.011;
aero.D1 = 0.01;
aero.D2 = 0.01;
rotor(2).aero = aero;

%Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
%CL_alpha
rotor(2).solidity = rotor(1).solidity; 
aero.cl_alpha = 5.156;         %1/rad - Lift slope, assumed NACA 0012 from graph
aero.Cd0 =  0.011;
aero.D1 = 0.0;
aero.D2 = 0.0;
rotor(2).aero = aero;


%% General params

params.kappaint = 1.28;
params.kappa = 1.15;

%% Testing

%Fcf= Prandtl_tip_loss(r,lambda,rotor);
%Prandtl_tip_loss(0.3,0,rotor(1))
%Prandtl_tip_loss([0,1],[2,3],rotor(1))

%lambda_bot = get_lambda_bot(F,r,pitch,rotor);
%lambda_bot = get_lambda_bot(1,0.3,0.2,rotor)
%lambda_bot = get_lambda_bot([0.58,0.4],[0,0.3],[0.2,0.2],rotor)

