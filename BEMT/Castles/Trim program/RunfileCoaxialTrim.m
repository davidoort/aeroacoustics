      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %    Non-linear, 6DOF  trim simulation for the Coaxial helicopter   %
      %                                                                   %
      %                     Master of science thesis                      %
      %                         February 2008                             %
      %                                                                   %
      %                     By John Campfens 1296434                      %
      %                 Delft University of Technology                    %
      %                 Faculty of Aerospace Engineering                  %
      %                                                                   %
      %                     Software is based on                          %
      %         Graduation project of Jasper van der Vorst (1998)         %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
warning off

% Define globals

global kts2ms deg2rad
% For Coaxialtrim
global V vpa h Vin

% Read data
Ka32data

% Initialisation and trim calculation
statesM=[];
controlsM=[];
coeffsM=[];

Vin=20;                                          % Initial forward speed
Vmax=80;                                        % Maximum forward speed [m/s]

for V = Vin:1:Vmax;                               % Defining airspeed range for trim

    vpa = 0*deg2rad;                            % Initial vertical path angle (up = pos) [rad]
    h   = 100*ft2m;                             % Initial height [m]

    [states,controls,coeffs] = CoaxialTrim;        % Initialise CoaxialTrim program
    coeffsM=[coeffsM; coeffs];                  % Putting coefficients into table
    statesM=[statesM; states'];                 % Putting states into table
    controlsM=[controlsM; controls];            % Putting controls into table

end

% Saving trim data
save coax_trimstates statesM controlsM coeffsM

% Initialize trim postprocessing
TrimOutput


