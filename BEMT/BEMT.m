%% Init

clear
close all
clc

addpath(genpath('./'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change parameters

coaxial.state.axial_vel = 0; %m/s 
coaxial.state.tangent_vel = 0; %m/s 
coaxial.state.trim = 1.0781;
coaxial.state.pitchdeg = 4;


epsilon = 0.0001; %convergence accuracy for Fcf and lambda



%% Run simple code

plots = true;
verbose = true;


if coaxial.state.tangent_vel == 0
    %BEMT_FF can also run with a tangential velocity of 0 but it will take
    %longer
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
    
    %[Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots); %for testing
    
else
    
    [Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots);
    
end

%% Iteration to trim the coaxial rotor and produce CT-CP plots

iter_pitchdeg = 0:1:18;

CT_arr = zeros(1,length(iter_pitchdeg));
CP_arr = zeros(1,length(iter_pitchdeg));

for idx = 1:length(iter_pitchdeg)
    tic
    [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,iter_pitchdeg(idx),'pitch_upper');
    toc
    if coaxial.state.CT < 0
        CT_arr(idx) = [];
        CP_arr(idx) = [];
        
    else
        CT_arr(idx) = coaxial.state.CT;
        CP_arr(idx) = coaxial.state.CP;
    end
    
end

plot(CP_arr,CT_arr)

%% Trim the coaxial rotor at a specified thrust coefficient

CT_desired = 0.005;

[collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT");


disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])
