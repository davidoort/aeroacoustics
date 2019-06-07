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
coaxial.state.pitchdeg = 14;


epsilon = 0.0001; %convergence accuracy for Fcf and lambda

% Testing

%[collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,coaxial.state.pitchdeg,'pitch_upper');


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

%% Iteration to trim the coaxial rotor and produce CT-CP validation plots

iter_pitchdeg = 0:1:18;

CT_arr = zeros(1,length(iter_pitchdeg));
CP_arr = zeros(1,length(iter_pitchdeg));
coaxial.state.axial_vel = 0; %m/s - comparison plots are for hover

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

% Verification plot H1 CP vs CT

if strcmpi(coaxial.name,"Harrington1")
    ax_xlim = [0 0.0006];
    ax_ylim = [0 0.007];
    if strcmpi(coaxial.type,"single")
        val_data =readmatrix('H1_single_fig2.csv');
    else
        val_data =readmatrix('H1_coax_fig2.csv');
    end
elseif strcmpi(coaxial.name,"Harrington2")
    ax_xlim = [0 0.001];
    ax_ylim = [0 0.01];
    if strcmpi(coaxial.type,"single")
        val_data =readmatrix('H2_single_fig3.csv');
    else
        val_data =readmatrix('H2_coax_fig3.csv');
    end
else
    val_data = [];
end
    
    
CP_exp = val_data(:,1);
CT_exp = val_data(:,2);

figure(1)
clf;
hold on
scatter(CP_exp,CT_exp)
plot(CP_arr,CT_arr)
xlabel('$C_P$','Interpreter','latex')
ylabel('$C_T$','Interpreter','latex')
xlim(ax_xlim)
ylim(ax_ylim)
title(strcat(coaxial.name, " ",coaxial.type," rotor"))
legend('Experiment','BEMT')

%% Verification plots with FVM for inflow, CT and CP distributions. Trim the coaxial rotor at a specified thrust coefficient

CT_desired = 0.004;

[collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT");


disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])


plots = true;
verbose = false;
[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);

%% Axial flight plots - not working yet, I need an efficient way of detecting stall


iter_pitchdeg_axial = [10]; %tweak this selection

axial_vel_range = 0:5:20;

CT_arr_axial = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));
CP_arr_axial = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));
axial_vel_arr = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));

coaxial.state.trim = 1; %possibly I can add a trim procedure later, but this at least ensures convergence fast
plots=false;
verbose=false;


for idx = 1:length(iter_pitchdeg)
    %coaxial.state.pitchdeg = iter_pitchdeg(idx);
    
    for axial_vel = axial_vel_range
        
        coaxial.state.axial_vel = axial_vel; %m/s 
        
        %[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        %coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
        
        tic
        [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,iter_pitchdeg(idx),'pitch_upper');
        toc
        
        if coaxial.state.CT < 0
            CT_arr_axial(idx) = [];
            CP_arr_axial(idx) = [];
            axial_vel_arr(idx) = [];
        else
            CT_arr_axial(idx) = coaxial.state.CT;
            CP_arr_axial(idx) = coaxial.state.CP;
            axial_vel_arr(idx) = axial_vel;
        end
    end
end


scatter(axial_vel_arr,CT_arr_axial)

