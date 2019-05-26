

addpath(genpath('./'))


atm = Atmosphere();
coaxial = Coaxial();


epsilon = 0.001; %convergence accuracy

plots = false;

if tangent_vel == 0
    %BEMT_FF can also run with a tangential velocity of 0 but it will take
    %longer
    [Thrust_ax, Torque_ax, Power_ax] = BEMT_axial(coaxial,atm,epsilon,plots);
    
    [Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots); %for testing
    
else
    
    [Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots);
    
end

