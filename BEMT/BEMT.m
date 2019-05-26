

addpath(genpath('./'))


atm = Atmosphere();
coaxial = Coaxial();


epsilon = 0.001; %convergence accuracy

plots = false;

[Thrust, Torque, Power] = BEMT_axial(coaxial,atm,epsilon,plots);

