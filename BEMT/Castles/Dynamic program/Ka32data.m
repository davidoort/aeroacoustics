% Ka32data provides all helicopter data of the Ka-32 coaxial helicopter. If
% necessary Puma data is used.

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

% Based on graduation project of Jasper van der Vorst (1998)

% Define globals
% For CoaxialTrim
global F0 W deg2rad
% For CoaxialModel
global Cla c R I_beta eps twist Kbeta
global F0 K_fus Vol_fus sigma
global alpha0_hs x_hs S_hs Cla_hs beta0_fin x_fin z_fin S_fin Cla_fin CDS k Cd
global n W Omega zhu yhu xhu m Iyy Ixx Jxz theta_tw
global Izz Nb tau_lambda xhl yhl zhl Vtip h_nd h0
% For Rho
global hstrat rho0 lambda T0 g R0

% Conversions
deg2rad = pi/180;       % Degrees to radians
kts2ms  = 1852/3600;    % Knots to meter per second
ft2m    = 0.3048;       % Feet to meter
hp2watt = 745.69987;    % HorsePower to watt

% Environmental constants
g	   = 9.80665; % Gravitational acceleration               [m/s^2]
R0	   = 287.05;  % Specific gas constant of air             [J/kg/K]
T0	   = 288.15;  % Sea level temperature in ISA             [K]
hstrat = 11000;   % Altitude at which stratosphere begins    [m]
rho0   = 1.2250;  % Sea level density in ISA                 [kg/m^3]
lambda = -0.0065; % Standard atmosphere temperature gradient [K/m]

% General helicopter parameters
m = 10000; % Helicopter mass   [kg]
W = m*g;  % Helicopter weight  [N]

Ixx = 9638;  % Helicopter moment of inertia about x-axis [kg*m^2]
Iyy = 33240; % Helicopter moment of inertia about y-axis [kg*m^2]
Izz = 25889; % Helicopter moment of inertia about z-axis [kg*m^2]
Jxz	= 2226;  % Heli product of inertia about x & z-axis  [kg*m^2]

dxcg = 0; % Displacement of cg along x-axis wrt ref point [m]
dycg = 0; % Displacement of cg along y-axis wrt ref point [m]
dzcg = 0; % Displacement of cg along z-axis wrt ref point [m]

n = 4.65; % (van Holten p51) [-]

% Rotor parameters (upper rotor and lower rotor)
Omega       = 28.4277;              % Angular velocity of rotor [rad/s]
R           = 7.95;                 % Main rotor radius [m]
Vtip        = 226;                  % Rotor tip speed [m/s]
c           = 0.48;                 % Rotor blade chord [m]
Nb          = 3;                    % Rotor number of blades [-]
Nr          = 2;                    % Number of rotors [-]
sigma       = Nb*c/(pi*R);          % Rotor solidity [-]
Cla         = 5.73;                 % Main rotor liftgradient of NACA 23012  [1/rad]
Cd          = 0.01;                 % Approximated drag coefficient
I_beta      = 1280;                 % Rotor blade flap mom of inertia [kg*m^2]
theta_tw    = -6*deg2rad;			% Rotor blade linear twist [rad]
e	        = 0;        			% Flapping hinge offset [m]
eps	        = e/R;          		% Rotor blade offset (e/R) [-]
tau_lambda  = 0.1;                  % Time constant                   [s]
Kbeta	    = 33032;                % Equivalent spring constant      [N*m/rad]
gamma       = rho0*Cla*c*R^4/I_beta;% Lock nr at sea level [-]
k           = 1.15;                 % (van Holten p30)                [-]
h_nd        = 0.189;                % Non-dimensional separation distance [-]
h0          = 0.189*R;              % Separation distance

% Determining positions of the rotors
xhu = 0;        % Hub x-position relative to cg [m]
yhu = 0;        % Hub y-position relative to cg [m]
zhu = 2.186+h0; % Hub z-position relative to cg [m]

xhl = 0;        % Hub x-position relative to cg [m]
yhl = 0;        % Hub y-position relative to cg [m]
zhl = 2.186;    % Hub z-position relative to cg [m]

% Fuselage
F0	    = 1.8;    % Parasite drag area of the helicopter    [m^2]
S_fus   = 26;     % Fuselage surface                        [m^2]
K_fus	= 0.83;	  % Correction coeff in fus pitching moment [-]
Vol_fus	= 6.11;   % Equivalent volume of circular body      [m^3]
CDS     = 4;      % Eq. flat plate area                     [m^2]
n       = 4.65;   % (van Holten p51)                        [-]

% Horizontal stabilizer
Cla_hs	  = 4;          % Horizontal stabilizer lift gradient       [1/rad]
alpha0_hs = -1.5*pi/180;% Horizontal stabilizer incidence           [rad]
S_hs	  = 1.335;      % Horizontal stabilizer area                [m^2]
x_hs	  = 8.92;       % Horizontal stabilizer x-positon rel to cg [m]
K_hs	  = 1.5;        % Horizontal stabilizer downwash factor     [-]
% 
% Vertical fin
Cla_fin	  = 4;         % Vertical fin lift gradient             [1/rad]
S_fin	  = 1.395;	   % Vertical fin area                      [m^2]
beta0_fin = -1*pi/180; % Vertical fin incidence                 [rad]
x_fin	  = 8.563;     % Vertical fin x-position relative to cg [m]
z_fin	  = 0;         % Vertical fin z-position relative to cg [m]
