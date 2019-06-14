%% Init

clear all
close all
clc

addpath(genpath('./'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change parameters

coaxial.state.axial_vel = 0; %m/s 
coaxial.state.tangent_vel = 40; %m/s 
coaxial.state.trim= 1;
coaxial.state.collective = 60; %collective in deg

epsilon = 0.0001; %convergence accuracy for Fcf and lambda -> 0.0001

warning('off')

% Testing 

plots= true;
verbose= true;
method='airfoil'; %'leishman','airfoil'

if strcmpi(method,'leishman')
    [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT_FF(coaxial,atm,epsilon,plots,verbose);
elseif strcmpi(method,'airfoil')
    [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT_iter(coaxial,atm,epsilon,plots,verbose);
end


%% Iteration to trim the coaxial rotor and produce CT-CP validation plots

iter_pitchdeg = 0:1:18;
method = 'leishman';
coaxial.state.axial_vel = 0; %m/s - comparison plots are for hover

for rotor_type = ["single","coaxial"] 
    
    CT_arr = zeros(1,length(iter_pitchdeg));
    CP_arr = zeros(1,length(iter_pitchdeg));
    %C_T_SMT_arr = zeros(1,length(iter_pitchdeg));
    C_P_SMT_arr = zeros(1,length(iter_pitchdeg));
    
    coaxial.type = rotor_type;

    for idx = 1:length(iter_pitchdeg)
        %if strcmpi(coaxial.type,'single')
            
        %elseif strcmpi(coaxial.type,'coaxial')
            
        %end
        tic
        [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,iter_pitchdeg(idx),'pitch_upper',method);
        toc
        if coaxial.state.CT < 0
            CT_arr(idx) = [];
            CP_arr(idx) = [];

        else
            CT_arr(idx) = coaxial.state.CT;
            CP_arr(idx) = coaxial.state.CP;
            
            if strcmpi(coaxial.type,"coaxial")
                % SMT Coaxial
                C_P_SMT_arr(idx) = 0.5*coaxial.params.kappaint*coaxial.params.kappa*(coaxial.state.CT)^(3/2)+(coaxial.rotor(2).solidity + coaxial.rotor(1).solidity)*coaxial.rotor(1).aero.Cd0/8; %using Cd0 of bottom rotor but could be for top rotor, this is for validation purposes
            else
                %SMT Single
                C_P_SMT_arr(idx) = coaxial.params.kappa*coaxial.state.CT^(3/2)/sqrt(2) + coaxial.rotor(1).solidity*coaxial.rotor(1).aero.Cd0/8;

            end
        end

    end

    % Verification plot H1 or H2 CP vs CT

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
    legend('-DynamicLegend');
    hold all
    scatter(CP_exp,CT_exp,'DisplayName',strcat("Experiment ", coaxial.type," rotor"))
    plot(CP_arr,CT_arr,'DisplayName',strcat("BEMT ",coaxial.type," rotor"),'LineWidth',2)
    plot(C_P_SMT_arr,CT_arr,'b-.','DisplayName',strcat("SMT ",coaxial.type," rotor"))
    xlabel('$C_P$','Interpreter','latex')
    ylabel('$C_T$','Interpreter','latex')
    xlim(ax_xlim)
    ylim(ax_ylim)
    title(string(coaxial.name))
    
end


%% Verification plots with FVM for inflow, CT and CP distributions. Trim the coaxial rotor at a specified thrust coefficient

CT_desired = 0.004;
method='leishman';


[collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT",method);


disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])


plots = true;
verbose = false;

%Don't change for now!
[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);

%% Axial flight plots - not working yet, I need an efficient way of detecting stall

%SOMETHING SUPER WEIRD IS GOING ON HERE


iter_pitchdeg_axial = [5]; %tweak this selection

axial_vel_range = 0:1:2;

CT_arr_axial = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));
CP_arr_axial = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));
axial_vel_arr = zeros(1,length(iter_pitchdeg_axial)*length(axial_vel_range));

coaxial.state.trim = 1; %possibly I can add a trim procedure later, but this at least ensures convergence fast
plots=false;
verbose=true;


for idx = 1:length(iter_pitchdeg_axial)
    %coaxial.state.pitchdeg = iter_pitchdeg(idx);
    
    for axial_vel = axial_vel_range
        
        coaxial.state.axial_vel = axial_vel; %m/s 
        
        %[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        %coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
        
        try
        tic
        [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,iter_pitchdeg_axial(idx),'pitch_upper');
        toc
        catch
            disp("Negative thrust")
            coaxial.state.CT = -1; %random
        end
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

