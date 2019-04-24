% M-file for outputting the yaw quickness of the coaxial helicopter 

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

close all
clear all
clc

% Loading quickness data (distrubution of Kupp=0.5 and Klow=-0.5)
load yawquickness10deg
load yawquickness20deg
load yawquickness30deg
load yawquickness40deg
load yawquickness50deg
load yawquickness60deg

% Loading quickness data (distrubution of Kupp=2.0 and Klow=-2.0) = sit. a
load yawquickness10dega
load yawquickness20dega
load yawquickness30dega
load yawquickness40dega
load yawquickness50dega
load yawquickness60dega

% Loading trim conditions
load coax_trimstates

% Defining intitial speed
Vin=20;

% Combining intial speed to tabulated data
N=Vin+1;

in_states=states(N,:)';
in_controls=controls(N,:)';

% Determining maximum yaw attitude and position of maximum
[ymax10,posmax10]=max(x10(:,7));
[ymax20,posmax20]=max(x20(:,7));
[ymax30,posmax30]=max(x30(:,7));
[ymax40,posmax40]=max(x40(:,7));
[ymax50,posmax50]=max(x50(:,7));
[ymax60,posmax60]=max(x60(:,7));

% Determining maximum yaw attitude and position of maximum (sit. a)
[ymax10a,posmax10a]=max(x10a(:,7));
[ymax20a,posmax20a]=max(x20a(:,7));
[ymax30a,posmax30a]=max(x30a(:,7));
[ymax40a,posmax40a]=max(x40a(:,7));
[ymax50a,posmax50a]=max(x50a(:,7));
[ymax60a,posmax60a]=max(x60a(:,7));

% Determining quicknesses
Q10=(max(x10(:,6))*180/pi)/(ymax10*180/pi)
Q20=(max(x20(:,6))*180/pi)/(ymax20*180/pi)
Q30=(max(x30(:,6))*180/pi)/(ymax30*180/pi)
Q40=(max(x40(:,6))*180/pi)/(ymax40*180/pi)
Q50=(max(x50(:,6))*180/pi)/(ymax50*180/pi)
Q60=(max(x60(:,6))*180/pi)/(ymax60*180/pi)

% Determining quicknesses (sit. a)
Q10a=(max(x10a(:,6))*180/pi)/(ymax10a*180/pi);
Q20a=(max(x20a(:,6))*180/pi)/(ymax20a*180/pi);
Q30a=(max(x30a(:,6))*180/pi)/(ymax30a*180/pi);
Q40a=(max(x40a(:,6))*180/pi)/(ymax40a*180/pi);
Q50a=(max(x50a(:,6))*180/pi)/(ymax50a*180/pi);
Q60a=(max(x60a(:,6))*180/pi)/(ymax60a*180/pi);

% Determining minimal yaw attitude
ymin10=(min(x10(posmax10:end,7)))*180/pi;
ymin20=(min(x20(posmax20:end,7)))*180/pi;
ymin30=(min(x30(posmax30:end,7)))*180/pi;
ymin40=(min(x40(posmax40:end,7)))*180/pi;
ymin50=(min(x50(posmax50:end,7)))*180/pi;
ymin60=(min(x60(posmax60:end,7)))*180/pi;

% Determining minimal yaw attitude (sit. a)
ymin10a=(min(x10a(posmax10a:end,7)))*180/pi;
ymin20a=(min(x20a(posmax20a:end,7)))*180/pi;
ymin30a=(min(x30a(posmax30a:end,7)))*180/pi;
ymin40a=(min(x40a(posmax40a:end,7)))*180/pi;
ymin50a=(min(x50a(posmax50a:end,7)))*180/pi;
ymin60a=(min(x60a(posmax60a:end,7)))*180/pi;


figure(1)
subplot(1,3,1)
plot(t10a,x10a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t10a,x10a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t10a,y10a(:,1)*180/pi,t10a,y10a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(2)
subplot(1,3,1)
plot(t20a,x20a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t20a,x20a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t20a,y20a(:,1)*180/pi,t20a,y20a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(3)
subplot(1,3,1)
plot(t30a,x30a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t30a,x30a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t30a,y30a(:,1)*180/pi,t30a,y30a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(4)
subplot(1,3,1)
plot(t40a,x40a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t40a,x40a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t40a,y40a(:,1)*180/pi,t40a,y40a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(5)
subplot(1,3,1)
plot(t50a,x50a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t50a,x50a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t50a,y50a(:,1)*180/pi,t50a,y50a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(6)
subplot(1,3,1)
plot(t60a,x60a(:,7)*180/pi)
legend('\psi_f')
xlabel('t [s]')
ylabel('fuselage yaw (heading) angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t60a,x60a(:,6)*180/pi)
legend('r')
xlabel('t [s]')
ylabel('yaw rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t60a,y60a(:,1)*180/pi,t60a,y60a(:,2)*180/pi)
legend('\theta_0_u','\theta_0_l')
xlabel('t [s]')
ylabel('Collective pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

figure(7)
plot(ymin10,Q10,'*',ymin20,Q20,'*',ymin30,Q30,'*',ymin40,Q40,'*',ymin50,Q50,'*',ymin60,Q60,'*')
axis([0 60 0 2.5])
legend('\psi_c_h_a_n_g_e = 10 deg','\psi_c_h_a_n_g_e = 20 deg','\psi_c_h_a_n_g_e = 30 deg','\psi_c_h_a_n_g_e = 40 deg','\psi_c_h_a_n_g_e = 50 deg','\psi_c_h_a_n_g_e = 60 deg')
xlabel('Minimum attitude change \Delta\psi_m_i_n [deg]')
ylabel('yaw attitude quickness [1/s]')
grid on

figure(8)
plot(ymin10a,Q10a,'*',ymin20a,Q20a,'*',ymin30a,Q30a,'*',ymin40a,Q40a,'*',ymin50a,Q50a,'*',ymin60a,Q60a,'*')
axis([0 61 0 2.5])
legend('\psi_c_h_a_n_g_e = 10 deg','\psi_c_h_a_n_g_e = 20 deg','\psi_c_h_a_n_g_e = 30 deg','\psi_c_h_a_n_g_e = 40 deg','\psi_c_h_a_n_g_e = 50 deg','\psi_c_h_a_n_g_e = 60 deg')
xlabel('Minimum attitude change \Delta\psi_m_i_n [deg]')
ylabel('yaw attitude quickness [1/s]')
grid on
