% M-file for outputting the dynamic response of the coaxial helicopter to a 
% disturbance of the collective pitch

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

clear all
close all
clc

load collectiveresponse

ft2m    = 0.3048;
m2ft    = 1/ft2m;
kts2ms  = 1852/3600;
deg2rad = pi/180;
rad2deg = 180/pi;

figure(1)
subplot(1,2,1)
plot(t,-x(:,12)*m2ft)
legend('h')
xlabel('t [s]')
ylabel('Altitude [ft]')
title('(a)')
grid on
subplot(1,2,2)
plot(t,-x(:,12)*m2ft)
legend('h')
xlabel('t [s]')
ylabel('Altitude [ft]')
xlim([0 15])
title('(b)')
grid on

figure(2)
subplot(1,3,1)
plot(t,x(:,8)*rad2deg)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
title('(a)')
subplot(1,3,2)
plot(t,x(:,9)*rad2deg)
legend('\phi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
title('(b)')
subplot(1,3,3)
plot(t,x(:,7)*rad2deg)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw angle [deg]')
grid on
title('(c)')


figure(3)
plot(t,-x(:,12)*1/ft2m)
legend('h')
xlabel('t [s]')
ylabel('altitude [ft]')
grid on

figure(4)
subplot(2,1,1)
plot(t,x(:,6)*rad2deg)
legend('p')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 15])
subplot(2,1,2)
plot(t,x(:,7)*rad2deg)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
xlim([0 15])

figure(5)
subplot(1,2,1)
plot(t,y(:,1)*rad2deg,t,y(:,2)*rad2deg,t,y(:,1)*rad2deg+y(:,2)*rad2deg)
legend('\theta_0_u','\theta_0_l','\theta_0_u+\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angles [deg]')
title('(a)')
%xlim([0 15])
grid on
subplot(1,2,2)
plot(t,y(:,3)*rad2deg,t,y(:,4)*rad2deg)
legend('\theta_1_s','\theta_1_c')
xlabel('t [s]')
ylabel('Cyclic pitch angles [deg]')
title('(b)')
grid on






