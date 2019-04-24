% M-file for outputting the roll quickness of the coaxial helicopter 

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

close all
clear all
clc

% Loading quickness data
load rollquickness10deg
load rollquickness20deg
load rollquickness30deg
load rollquickness40deg


% Loading trim conditions
load coax_trimstates

% Defining intitial speed
Vin=20;

% Combining intial speed to tabulated data
N=Vin+1;

in_states=states(N,:)';
in_controls=controls(N,:)';

% Determining maximum roll attitude and position of maximum
[pmax10,posmax10]=max(x10(:,9));
[pmax20,posmax20]=max(x20(:,9));
[pmax30,posmax30]=max(x30(:,9));
[pmax40,posmax40]=max(x40(:,9));
% [pmax50,posmax50]=max(x50(:,9));
% [pmax60,posmax60]=max(x60(:,9));

% Determining quicknesses
Q10=(max(x10(:,4))*180/pi)/(pmax10*180/pi);
Q20=(max(x20(:,4))*180/pi)/(pmax20*180/pi);
Q30=(max(x30(:,4))*180/pi)/(pmax30*180/pi);
Q40=(max(x40(:,4))*180/pi)/(pmax40*180/pi);
% Q50=(max(x50(:,4))*180/pi)/(pmax50*180/pi);
% Q60=(max(x60(:,4))*180/pi)/(pmax60*180/pi);

% Determining minimal roll attitude
pmin10=(min(x10(posmax10:end,9)))*180/pi;
pmin20=(min(x20(posmax20:end,9)))*180/pi;
pmin30=(min(x30(posmax30:end,9)))*180/pi;
pmin40=(min(x40(posmax40:end,9)))*180/pi;
% pmin50=(min(x50(posmax50:end,9)))*180/pi;
% pmin60=(min(x60(posmax60:end,9)))*180/pi;

figure(1)
subplot(1,3,1)
plot(t10,x10(:,9)*180/pi)
legend('\phi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t10,x10(:,4)*180/pi)
legend('p')
xlabel('t [s]')
ylabel('roll rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t10,y10(:,4)*180/pi)
legend('\theta_1_c')
xlabel('t [s]')
ylabel('Lateral cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')


figure(2)
subplot(1,3,1)
plot(t20,x20(:,9)*180/pi)
legend('\phi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t20,x20(:,4)*180/pi)
legend('p')
xlabel('t [s]')
ylabel('roll rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t20,y20(:,4)*180/pi)
legend('\theta_1_c')
xlabel('t [s]')
ylabel('Lateral cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(3)
subplot(1,3,1)
plot(t30,x30(:,9)*180/pi)
legend('\phi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t30,x30(:,4)*180/pi)
legend('p')
xlabel('t [s]')
ylabel('roll rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t30,y30(:,4)*180/pi)
legend('\theta_1_c')
xlabel('t [s]')
ylabel('Lateral cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(4)
subplot(1,3,1)
plot(t40,x40(:,9)*180/pi)
legend('\phi_f')
xlabel('t [s]')
ylabel('fuselage roll angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t40,x40(:,4)*180/pi)
legend('p')
xlabel('t [s]')
ylabel('roll rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t40,y40(:,4)*180/pi)
legend('\theta_1_c')
xlabel('t [s]')
ylabel('Lateral cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(5)
plot(pmin10,Q10,'*',pmin20,Q20,'*',pmin30,Q30,'*',pmin40,Q40,'*')
axis([0 50 0 2.5])
legend('\phi_c_h_a_n_g_e = 10 deg','\phi_c_h_a_n_g_e = 20 deg','\phi_c_h_a_n_g_e = 30 deg','\phi_c_h_a_n_g_e = 40 deg')
xlabel('Minimum attitude change \Delta\phi_m_i_n [deg]')
ylabel('roll attitude quickness [1/s]')
grid on


