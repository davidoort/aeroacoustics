
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
coaxial.state.trim = 1;
coaxial.state.pitchdeg = 8.71;


epsilon = 0.001; %convergence accuracy



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

k = 0.001;
eps = 0.0001;
plots = false;
verbose = false;

iter_pitchdeg = 1:0.1:18;

CT_arr = zeros(1,length(iter_pitchdeg));
CP_arr = zeros(1,length(iter_pitchdeg));

for i = 1:length(iter_pitchdeg)
    
    coaxial.state.pitchdeg = iter_pitchdeg(i);
    
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
    
    while abs(coaxial.state.net_torque)>eps
        coaxial.state.trim = coaxial.state.trim + k*coaxial.state.net_torque;
        
        [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
    
    end
    
    if coaxial.state.CT < 0
        CT_arr(i) = [];
        CP_arr(i) = [];
        
    else
        CT_arr(i) = coaxial.state.CT;
        CP_arr(i) = coaxial.state.CP;
    end
    
    
    
end

plot(CP_arr,CT_arr)


%% Iterate to trim the coaxial rotor at a specified thrust coefficient

k1 = 50;
k2 = 0.001;
eps1 = 0.00001;
eps2 = 0.001;

plots = false;
verbose = false;

CT = 0.004;


[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);

thrust_error = CT - coaxial.state.CT; 

while abs(thrust_error) > eps1
    
    coaxial.state.pitchdeg = coaxial.state.pitchdeg + k1*thrust_error;
    
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
    
    while abs(coaxial.state.net_torque)>eps2
        coaxial.state.trim = coaxial.state.trim + k2*coaxial.state.net_torque;
        
        [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
    
    end
    
    thrust_error = CT - coaxial.state.CT; 
    
end

