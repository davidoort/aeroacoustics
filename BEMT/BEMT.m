

addpath(genpath('./'))


atm = Atmosphere();
coaxial = Coaxial();


epsilon = 0.001; %convergence accuracy

plots = true;

[Thrust, Torque, Power] = BEMT_axial(coaxial,atm,epsilon,plots);

