% CoaxialTrim determines the trim conditions of the coaxial helicopter
% following the Newton Iteration method.

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

% Based on graduation project of Jasper van der Vorst (1998)


function [states,controls,coeffs] = CoaxialTrim

% Define globals
global V vpa W h deg2rad R Omega F0 Vin
% For Coaxialmodel
global x01 Vtip

if V==Vin
% Initial guesses
D	       = 1/2*Rho(h)*V^2*F0;     % Drag                           [N]
theta	   = asin(-D*cos(vpa)/W);   % Helicopter pitch angle         [rad]
psi	       = 0;                     % Helicopter yaw angle           [rad]
phi	       = -0.008;                % Helicopter roll angle          [rad]
u	       = V*cos(theta);          % Airspeed along x-axis          [m/s]
v	       = V*sin(theta)*sin(phi); % Airspeed along y-axis          [m/s]
w	       = V*sin(theta)*cos(phi); % Airspeed along z-axis          [m/s]
p	       = 0;                     % Roll rate                      [rad/s]
q	       = 0;                     % Pitch rate                     [rad/s]
r	       = 0;                     % Yaw rate                       [rad/s]
x          = 0;                     % Position along Earth x-axis    [m]
y	       = 0;                     % Position along Earth y-axis    [m]
z	       = -h;                    % Position along Earth z-axis    [m]
L0_u       = 0.0453;                  % Norm unif induc downwash of upper rotor [-]
L0_l       = 0.0231;                  % Norm unif induc downwash of lower rotor [-]
theta_1c   = 0.05;                     % Lateral cyclic position        [deg]
theta_1s   = 0;                     % Longitudinal cyclic position   [deg]
theta_0u   = 0.22;                  % Upper rotor collective position [deg]
theta_0l   = 0.19;                  % Lower rotor collective position [deg]
rhoh       = Rho(h);

else    % Taking previous trim condition as initial condition for the next iteration step

D          = 1/2*Rho(h)*V^2*F0;     % Drag                           [N]
theta	   = x01(8);                % Helicopter pitch angle         [rad]
psi	       = 0;                     % Helicopter yaw angle           [rad]
phi	       = 0;                     % Helicopter roll angle          [rad]
u	       = V*cos(theta);          % Airspeed along x-axis          [m/s]
v	       = V*sin(theta)*sin(phi); % Airspeed along y-axis          [m/s]
w	       = V*sin(theta)*cos(phi); % Airspeed along z-axis          [m/s]
p	       = 0;                     % Roll rate                      [rad/s]
q	       = 0;                     % Pitch rate                     [rad/s]
r	       = 0;                     % Yaw rate                       [rad/s]
x          = 0;                     % Position along Earth x-axis    [m]
y	       = 0;                     % Position along Earth y-axis    [m]
z	       = -h;                    % Position along Earth z-axis    [m]
L0_u       = x01(13);               % Norm unif induc downwash of upper rotor [-]
L0_l       = x01(14);               % Norm unif induc downwash of lower rotor [-]
theta_1c   = 0;                     % Lateral cyclic position        [deg]
theta_1s   = x01(17);               % Longitudinal cyclic position   [deg]
theta_0u   = x01(15);               % Upper rotor collective position         [deg]
theta_0l   = x01(16);               % Lower rotor collective position         [deg]
rhoh	   = Rho(h);                % Density of air                 [kg/m^3]

end

states = [u v w p q r psi theta phi x y z L0_u L0_l]'; % Making state vector

trimvar	= [states(8)        % 1:  theta
           states(9)        % 2:  phi
		   theta_0u         % 3:  Collective of upper rotor
		   theta_0l         % 4:  Collective of lower rotor
		   theta_1s         % 5: longitudinal cyclic
		   theta_1c         % 6: lateral cyclic
           states(13)       % 7: lambda upper rotor
           states(14)];     % 8: Lambda lower rotor
		   
f  = 1;
nn = 0;

delta = 1e-11*(ones(length(trimvar),1));
dfdx  = zeros(length(trimvar));

% Initiating trim process

while max(abs(f)) > 1e-8
    nn 	       = nn+1;
        
    states(8)  = trimvar(1);
    states(9)  = trimvar(2);
    theta_0u   = trimvar(3);
    theta_0l   = trimvar(4);
    theta_1s   = trimvar(5);
    theta_1c   = trimvar(6);
    states(13) = trimvar(7);
    states(14) = trimvar(8);

    dot = CoaxialModel(states,theta_0u,theta_0l,theta_1s,theta_1c);

    f	  = [dot(1:6) 
             dot(13:14)];  % Making state vector f which has to equal zero

    oldtrimvar = trimvar;

	for i = 1:length(trimvar)
		perturb	   = zeros(length(trimvar),1);
		perturb(i) = delta(i);
		trimvar	   = oldtrimvar+perturb;

		states(8)  = trimvar(1);
        states(9)  = trimvar(2);
        theta_0u   = trimvar(3);
        theta_0l   = trimvar(4);
        theta_1s   = trimvar(5);
        theta_1c   = trimvar(6);
        states(13) = trimvar(7);
        states(14) = trimvar(8);

		dot = CoaxialModel(states,theta_0u,theta_0l,theta_1s,theta_1c);

% Executing Newton Iteration method        
        fnew     = [dot(1:6)
                    dot(13:14)];

  		df	      = (fnew-f)/delta(i);  % Calculating particial derivative
  		dfdx(:,i) = df;                 % Making Jacobian Matrix
       
    end

    inc	    = -inv(dfdx)*f;                      
    trimvar	= inc+trimvar;              % Determining new trim condition
end

% New trim variables in trimmed condition
    states(8)  = trimvar(1); 
    states(9)  = trimvar(2);
    theta_0u   = trimvar(3);
    theta_0l   = trimvar(4);
    theta_1s   = trimvar(5);
    theta_1c   = trimvar(6);
    states(13) = trimvar(7);
    states(14) = trimvar(8);

% Outputting trim variables
disp('Number of iterations:')
disp([nn])
disp('Solution:')
disp(['u [m/s]'])
disp([states(1)])
disp(['v [m/s]'])
disp([states(2)])
disp(['w [m/s]'])
disp([states(3)])
disp(['theta_f [deg]'])
disp([states(8)/deg2rad])
disp(['phi_f [deg]'])
disp([states(9)/deg2rad])
disp(['Li_up [-]'])
disp([states(13)])
disp(['Li_lo [-]'])
disp([states(14)])
disp([' theta_0u   theta_0l   theta_1s    theta_1c'])
disp([theta_0u  theta_0l        theta_1s       theta_1c])

f = CoaxialModel(states,theta_0u,theta_0l,theta_1s,theta_1c);

dots = f';

coeffs=[f(15) f(16) f(17) f(18)];

controls = [theta_0u theta_0l theta_1s theta_1c];

x01 = [states' controls];
