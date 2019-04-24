% TrimOutput determines the wanted outputs of the trim conditions

% Made by John Campfens (1296434)
% Delft University of Technology, February 2008

clear all
clc
close all

% Loading tabelized trimstates and controls
load coax_trimstates       % Trim states of coaxial helicopter model
load coax_trimstates_ni    % Trim states of coaxial helicopter model without interference 

% Loading Puma trim data
load trimarray              
load nlr_theta0_r.dat
load nlr_lat_r.dat
load nlr_lon_r.dat
load nlr_pitch_r.dat
load nlr_roll_r.dat

Ka32data

W2HP=0.001341;
kt2m = 0.5151;

E=0:1:80;   % Forward speed vector

% Plotting the different trimmed states and controls
figure(1);
plot(E,controls(:,1).*180./pi,E,controls(:,2).*180./pi,E,controlsM(:,1).*180./pi,M_V,M_FCS(:,1),'m-^',nlr_theta0_r(:,1).*kt2m,nlr_theta0_r(:,2),'k*');
legend('Upper rotor','Lower rotor','No interference','MoGeHM Puma trim','SAAF Puma measurements');
xlabel('Forward speed [m/s]')
ylabel('\theta_0 [deg]')
grid on
axis([0 70 6 16]);

figure(2);
plot(E,controls(:,4).*180./pi,E,controlsM(:,4).*180./pi,M_V,M_FCS(:,2),'m-^',nlr_lat_r(:,1).*kt2m,nlr_lat_r(:,2),'k*');
legend('Coaxial rotor','No interference','MoGeHM Puma trim','SAAF Puma measurements');
xlabel('Forward speed [m/s]')
ylabel('\theta_1_c [deg]')
grid on
axis([0 70 -0.5 3]);

figure(3);
plot(E,controls(:,3).*180./pi,E,controlsM(:,3).*180./pi,M_V,-M_FCS(:,3)+2.5,'m-^',nlr_lon_r(:,1).*kt2m,-nlr_lon_r(:,2)+5,'k*');
legend('Coaxial rotor','No interference','MoGeHM Puma trim','SAAF Puma measurements');
xlabel('Forward speed [m/s]')
ylabel('\theta_1_s [deg]')
grid on
axis([0 70 0 8]);

figure(4);
plot(E,states(:,8).*180./pi,E,statesM(:,8).*180./pi,M_V, M_u(:,1)*180/pi-2.4,'m-^',nlr_pitch_r(:,1).*kt2m,nlr_pitch_r(:,2)-2.4,'k*');
legend('Coaxial rotor','No interference','MoGeHM Puma trim','SAAF Puma measurements');
xlabel('Forward speed [m/s]')
ylabel('\theta_f [deg]')
grid on
axis([0 70 -6 1]);

figure(5);
plot(E,states(:,9).*180./pi,E,statesM(:,9).*180./pi,M_V, -M_u(:,2)*180/pi,'m-^',nlr_roll_r(:,1).*kt2m,nlr_roll_r(:,2),'k*');
legend('Coaxial rotor','No interference','MoGeHM Puma trim','SAAF Puma measurements');
xlabel('Forward speed [m/s]')
ylabel('\phi_f [deg]')
grid on
axis([0 70 -0.5 3]);

figure(6);
plot(E,states(:,13).*Vtip,E,states(:,14).*Vtip,E,statesM(:,13).*Vtip);%.*2,E,(states(:,13)+states(:,14)).*Vtip)
legend('Upper rotor','Lower rotor','No interference');
xlabel('Forward speed [m/s]')
ylabel('v_i [m/s]')
grid on
axis([0 70 0 12]);


figure(8);
subplot(2,1,1);
chi=atan((states(:,1)./226)./(states(:,13)-(states(:,3)./226)));
plot(E,chi.*180/pi)
legend('\chi');
xlabel('Forward speed [m./s]');
ylabel('Wake skew angle [deg]');
grid on
axis([0 70 0 90]);
subplot(2,1,2);
plot(E,controls(:,1).*180./pi,E,controls(:,2).*180./pi)
legend('Upper rotor','Lower rotor');
xlabel('Forward speed [m/s]');
ylabel('\theta_0 [deg]');
grid on
axis([0 70 6 16]);



%---------------------- Determining power curves ------------------------%

% Needed data for determining power curves


rhoh=1.225;
gamma = rhoh*Cla*c*R^2/I_beta;
sigma=Nb*c/(pi*R);
q=0;
dimless = rhoh*Vtip^2*pi*R^2;
E=0:1:80;
mu=(E/Vtip)';
alpha_cp = -atan2(states(:,3),states(:,1))+controls(:,3);
Cdp=0.02;

% Flapping and thrust coefficients needed for power calculations
wb=sqrt(Kbeta/I_beta)/Omega;
vb=sqrt(1+wb^2);

a0_u=(gamma./8*(1+(states(:,1)./Vtip).^2-2.*states(:,6)./Omega).*controls(:,1)-gamma./6.*(states(:,2)./Vtip).*controls(:,4)-gamma./6.*(states(:,1)./Vtip).*controls(:,3)+...
    gamma./8.*(4./5-8./5.*states(:,6)./Omega+2./3.*(states(:,1)./Vtip).^2).*theta_tw-gamma./6.*(-(states(:,3)./Vtip)+states(:,13))+gamma./12.*(states(:,1)./Vtip).*states(:,4)./Omega+...
    gamma./12.*(states(:,2)./Vtip).*states(:,5)./Omega)./(vb.^2-2.*states(:,6)./Omega);
a1_u=(gamma./6.*(states(:,2)./Vtip).*a0_u+gamma./3.*(states(:,1)./Vtip).*controls(:,1)-gamma./8.*(states(:,1)./Vtip).*(states(:,2)./Vtip).*controls(:,4)-...
    gamma./8.*(1-2.*states(:,6)./Omega+3./2.*(states(:,1)./Vtip).^2).*controls(:,3)+gamma./4.*(states(:,1)./Vtip).*theta_tw-gamma./4.*(states(:,1)./Vtip).*(-(states(:,3)./Vtip)+states(:,13))+...
    gamma./8.*states(:,4)./Omega-2.*states(:,5)./Omega)./(1-states(:,6)./Omega-(states(:,1)./Vtip).^2./2+(states(:,2)./Vtip).^2./2);

Ct_gl_u = 2.*states(:,13).*sqrt((mu.*cos(alpha_cp-a1_u)).^2+(mu.*sin(alpha_cp-a1_u)+states(:,13)).^2);

a0_l=(gamma./8.*(1+(states(:,1)./Vtip).^2-2.*states(:,6)./Omega).*controls(:,2)-gamma./6.*(states(:,2)./Vtip).*controls(:,4)-gamma./6.*(states(:,1)./Vtip).*controls(:,3)+...
    gamma./8.*(4./5-8./5.*states(:,6)./Omega+2./3.*(states(:,1)./Vtip).^2).*theta_tw-gamma./6.*(-(states(:,3)./Vtip)+states(:,14))+gamma./12.*(states(:,1)./Vtip).*states(:,4)./Omega+...
    gamma./12.*(states(:,2)./Vtip).*states(:,5)./Omega)./(vb.^2-2.*states(:,6)./Omega);
a1_l=(gamma./6.*(states(:,2)./Vtip).*a0_u+gamma./3.*(states(:,1)./Vtip).*controls(:,2)-gamma./8.*(states(:,1)./Vtip).*(states(:,2)./Vtip).*controls(:,4)-...
    gamma./8.*(1-2.*states(:,6)./Omega+3./2.*(states(:,1)./Vtip).^2).*controls(:,3)+gamma./4.*(states(:,1)./Vtip).*theta_tw-gamma./4.*(states(:,1)./Vtip).*(-(states(:,3)./Vtip)+states(:,14))+...
    gamma./8.*states(:,4)./Omega-2.*states(:,5)./Omega)./(1-states(:,6)./Omega-(states(:,1)./Vtip).^2./2+(states(:,2)./Vtip).^2./2);

Ct_gl_l = 2.*states(:,14).*sqrt((mu.*cos(alpha_cp-a1_l)).^2+(mu.*sin(alpha_cp-a1_l)+states(:,14)).^2);

a0_ni=(gamma./8.*(1+(statesM(:,1)./Vtip).^2-2.*statesM(:,6)./Omega).*controlsM(:,2)-gamma./6.*(statesM(:,2)./Vtip).*controlsM(:,4)-gamma./6.*(statesM(:,1)./Vtip).*controlsM(:,3)+...
    gamma./8.*(4./5-8./5.*statesM(:,6)./Omega+2./3.*(statesM(:,1)./Vtip).^2).*theta_tw-gamma./6.*(-(statesM(:,3)./Vtip)+statesM(:,14))+gamma./12.*(statesM(:,1)./Vtip).*statesM(:,4)./Omega+...
    gamma./12.*(statesM(:,2)./Vtip).*statesM(:,5)./Omega)./(vb.^2-2.*statesM(:,6)./Omega);
a1_ni=(gamma./6.*(statesM(:,2)./Vtip).*a0_u+gamma./3.*(statesM(:,1)./Vtip).*controlsM(:,2)-gamma./8.*(statesM(:,1)./Vtip).*(statesM(:,2)./Vtip).*controlsM(:,4)-...
    gamma./8.*(1-2.*statesM(:,6)./Omega+3./2.*(statesM(:,1)./Vtip).^2).*controlsM(:,3)+gamma./4.*(statesM(:,1)./Vtip).*theta_tw-gamma./4.*(statesM(:,1)./Vtip).*(-(statesM(:,3)./Vtip)+statesM(:,14))+...
    gamma./8.*statesM(:,4)./Omega-2.*statesM(:,5)./Omega)./(1-statesM(:,6)./Omega-(statesM(:,1)./Vtip).^2./2+(statesM(:,2)./Vtip).^2./2);

Ct_gl_ni = 2.*statesM(:,14).*sqrt((mu.*cos(alpha_cp-a1_ni)).^2+(mu.*sin(alpha_cp-a1_ni)+statesM(:,14)).^2);

% Parasite power
Ps=2*CDS.*.5.*rhoh.*E.^3;
Ps=Ps';

% Profile drag power
PpPd = 2./8.*sigma.*Cdp.*rhoh.*Vtip.^3.*pi.*R.^2.*(1+n.*(states(:,1)./Vtip).^2);
PpPd_ni= 2./8.*sigma.*Cdp.*rhoh.*Vtip.^3.*pi.*R.^2.*(1+n.*(statesM(:,1)./Vtip).^2);

% Induced power [W]
Pi_u = abs(k.*(Ct_gl_u.*dimless).*(states(:,13).*Vtip));
Pi_l = abs(k.*(Ct_gl_l.*dimless).*(states(:,14).*Vtip));
Pi=Pi_u+Pi_l;
Pi_ni=abs(k.*(Ct_gl_ni.*dimless).*(statesM(:,13).*Vtip));

Preq = Ps+Pi+PpPd;
Preq_ni=Ps+PpPd_ni+Pi_ni;

PavM=ones(81).*4380.*0.95.*0.87;
Pav=PavM(:,1);

figure(9)
plot(E.*3.6,Ps.*W2HP,E.*3.6,PpPd.*W2HP,E.*3.6,Pi.*W2HP,E.*3.6,Preq.*W2HP,E.*3.6,Pav)
xlabel('Forward speed [km./h]')
ylabel('Power [hp]')
legend('P_p_a_r','P_pP_d','P_i','P_r_e_q','P_a_v')
grid on
axis([0 250 0 4000]); 

figure(10);
plot(E.*3.6,Preq.*W2HP,E.*3.6,Preq_ni.*W2HP,E.*3.6,Pav)
xlabel('Forward speed [km./h]')
ylabel('Power [hp]')
legend('P_r_e_q','P_r_e_q_n_i','P_a_v')
grid on
axis([0 250 0 4000]);

