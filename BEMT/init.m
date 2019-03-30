
%% Atmospheric parameters

atm.rho = 1.225; %kg/m^3 sea-level air density
atm.gamma_air = 1.4; %[-] for air (diatomic gas value)
atm.R_air = 8.314510/0.0289645; %J/kg/K= N m / N s^2 m^-1 /K = m^2/s^2/K
atm.Temp = 20; %?C 
atm.c0 = sqrt(atm.gamma_air*atm.R_air*(atm.Temp+273.15)); %m/s - speed of sound

%% Coaxial rotor parameters
% Top rotor
rotor(1).sig_u = 0.0746;   % Blade stiffness
rotor(1).a_u = 2*pi;         % Lift coefficient
rotor(1).Nb = 2;             % Number of blades

% Bottom rotor
rotor(2).sig_l = rotor(1).sig_u; 
rotor(2).a_l = rotor(1).a_u;
rotor(2).rd = 0.82;  % annulus
rotor(2).A_Ad = 1/rotor(2).rd^2;

% Drag coefficients
Cd0 = 0.011;
D1 = 0.01;
D2 = 0.01;