

addpath(genpath('./'))

%% Instantiate objects

atm = Atmosphere();
coaxial = Coaxial();

%% Change parameters

coaxial.state.airspeed = 5; %m/s 
coaxial.state.trim = 1.4;
%coaxial.state.incidence_deg = 20; %deg - positive downward

coaxial.updateDependenciesPU

epsilon = 0.001; %convergence accuracy

plots = true;

%% Run code

if coaxial.state.tangent_vel == 0
    %BEMT_FF can also run with a tangential velocity of 0 but it will take
    %longer
    [Thrust_ax, Torque_ax, Power_ax] = BEMT_axial(coaxial,atm,epsilon,plots);
    
    %[Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots); %for testing
    
else
    
    [Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots);
    
end

