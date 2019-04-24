% M-file for outputing the maximum pitch attitude quickness
close all
clear all
clc

% Loading quickness data
load maxpitchquickness5deg
load maxpitchquickness10deg
load maxpitchquickness15deg
load maxpitchquickness20deg
load maxpitchquickness25deg
load maxpitchquickness30deg

% Loading trim conditions
load coax_trimstates

% Defining intitial speed
Vin=20;

% Combining intial speed to tabulated data
N=Vin+1;

in_states=states(N,:)';
in_controls=controls(N,:)';

% Determining maximum pitch attitude and position of maximum
[pmax5,posmax5]=max(x5(:,8));
[pmax10,posmax10]=max(x10(:,8));
[pmax15,posmax15]=max(x15(:,8));
[pmax20,posmax20]=max(x20(:,8));
[pmax25,posmax25]=max(x25(:,8));
[pmax30,posmax30]=max(x30(:,8));

% Determining quicknesses
Q5=(max(x5(:,5))*180/pi)/((pmax5-in_states(8))*180/pi);
Q10=(max(x10(:,5))*180/pi)/((pmax10-in_states(8))*180/pi);
Q15=(max(x15(:,5))*180/pi)/((pmax15-in_states(8))*180/pi);
Q20=(max(x20(:,5))*180/pi)/((pmax20-in_states(8))*180/pi);
Q25=(max(x25(:,5))*180/pi)/((pmax25-in_states(8))*180/pi);
Q30=(max(x30(:,5))*180/pi)/((pmax30-in_states(8))*180/pi);

% Determining minimal pitch attitude
pmin5=(min(x5(posmax5:end,8))-in_states(8))*180/pi;
pmin10=(min(x10(posmax10:end,8))-in_states(8))*180/pi;
pmin15=(min(x15(posmax15:end,8))-in_states(8))*180/pi;
pmin20=(min(x20(posmax20:end,8))-in_states(8))*180/pi;
pmin25=(min(x25(posmax25:end,8))-in_states(8))*180/pi;
pmin30=(min(x30(posmax30:end,8))-in_states(8))*180/pi;

figure(1)
subplot(1,3,1)
plot(t5,(x5(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t5,x5(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t5,y5(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')
% 
% % figure(2)
% % plot(t5,y5(:,3)*180/pi)
% 
% 
figure(2)
subplot(1,3,1)
plot(t10,(x10(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t10,x10(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t10,y10(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')
% 
figure(3)
subplot(1,3,1)
plot(t15,(x15(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t15,x15(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t15,y15(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')
% 
figure(4)
subplot(1,3,1)
plot(t20,(x20(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t20,x20(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t20,y20(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')
% 
figure(5)
subplot(1,3,1)
plot(t25,(x25(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t25,x25(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t25,y25(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')
% 
figure(6)
subplot(1,3,1)
plot(t30,(x30(:,8)-in_states(8))*180/pi)
legend('\theta_f')
xlabel('t [s]')
ylabel('fuselage pitch angle [deg]')
grid on
xlim([0 6])
title('(a)')
subplot(1,3,2)
plot(t30,x30(:,5)*180/pi)
legend('q')
xlabel('t [s]')
ylabel('pitch rate [deg/s]')
grid on
xlim([0 6])
title('(b)')
subplot(1,3,3)
plot(t30,y30(:,3)*180/pi)
legend('\theta_1_s')
xlabel('t [s]')
ylabel('Longitudinal cyclic pitch angle [deg]')
grid on
xlim([0 6])
title('(c)')

% figure(2)
% plot(pmin5,Q5,'*',pmin10,Q10,'*',pmin15,Q15,'*',pmin20,Q20,'*',pmin25,Q25,'*',pmin30,Q30,'*')
% axis([0 30 0 3])
% legend('\theta_c_h_a_n_g_e = 5 deg','\theta_c_h_a_n_g_e = 10 deg','\theta_c_h_a_n_g_e = 15 deg','\theta_c_h_a_n_g_e = 20 deg','\theta_c_h_a_n_g_e = 25 deg','\theta_c_h_a_n_g_e = 30 deg')
% xlabel('Minimum attitude change \Delta\theta_m_i_n [deg]')
% ylabel('pitch attitude quickness [1/s]')
% grid on
% 
% figure(3)
% plot(t25,y25(:,3)*180/pi,t30,y30(:,3)*180/pi)
% % figure(2)
% % plot(t5,x5(:,5))
% % 
% % 
% % figure(3)
% % plot(t5,y5(:,3))
