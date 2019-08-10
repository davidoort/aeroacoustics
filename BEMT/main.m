%% Init

clear
close all
clc

addpath(genpath('./'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change parameters

coaxial.state.axial_vel = 10; %m/s 
coaxial.state.forward_vel = 10; %m/s 
coaxial.state.side_vel = 10; %m/s
coaxial.state.trim = 1; %collective_l/collective_u
coaxial.state.collective = 20; %collective in deg - geometric pitch angle at the root of the UPPER rotor blades!

epsilon = 0.0001; %convergence accuracy for Prandtl tip function and inflow ratio

%warning('off')

% Testing 

plots= true;
verbose= true;
debug = false;
method='leishman'; %'leishman','airfoil'

%tic
%[Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
%toc

[Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);

%legend below
%{
%Power is a 1x2 matrix -> [P_upper; P_lower] [W]

%Forces is a 3x2 matrix -> [Fx_upper, Fy_upper, Fz_upper; Fx_lower, Fy_lower, Fz_lower] [N]
%using the right-handed coordinate system with x pointing forward, z pointing up

%Moments is a 3x2 matrix -> [Mx_upper, My_upper, Mz_upper; Mx_lower, My_lower, Mz_lower] [Nm] 
%using the same coordinate system as Forces

%CT is the thrust coefficient 1x2 matrix -> [CT_u; CT_l];
%CP is the torque coefficient 1x2 matrix -> [CP_u;CP_l];

%} 

%%%%% GENERATE DIFFERENCE (between methods) DISK PLOTS SOMEWHERE!!!!! ->
%%%%% this would be nice in the article (especially if I can explain where the biggest differences come from)
%% Iteration to trim the coaxial rotor and produce CT-CP validation plots

iter_pitchdeg = 0:1:18;
method = 'leishman'; %airfoil took about 6 mins to run

coaxial.state.forward_vel = 0; %m/s  - comparison plots are for hover
SMT = false;
%for method = ["leishman","single"]

for axial_vel = 0
    coaxial.state.axial_vel = axial_vel; %m/s - comparison plots are for hover
    rotors = ["single","coaxial"];
    for rotor_type = rotors
        %%
        
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
        set(gca,'FontSize',16)
        scatter(CP_exp,CT_exp,'DisplayName',strcat("Experiment ", coaxial.type," rotor"))
        plot(CP_arr,CT_arr,'DisplayName',strcat("BEMT ",coaxial.type," rotor"),'LineWidth',2)
        if SMT
            plot(C_P_SMT_arr,CT_arr,'b-.','DisplayName',strcat("SMT ",coaxial.type," rotor"))
        end
        xlabel('$C_P$','Interpreter','latex')
        ylabel('$C_T$','Interpreter','latex')
        xlim(ax_xlim)
        ylim(ax_ylim)
        title(string(coaxial.name))
        
    end
end
%% Verification plots with FVM for inflow, CT and CP distributions. Trim the coaxial rotor at a specified thrust coefficient

CT_desired = 0.004;
method='leishman';

coaxial.state.axial_vel = 0; %m/s - hover
coaxial.state.forward_vel = 0; %m/s - hover
epsilon = 0.001; %convergence accuracy for Fcf and lambda -> 0.0001

[collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT",method);


disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])


plots = true;
verbose = false;
debug = false;
%Don't change for now!
[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
%% Axial flight plots - not working yet, I need an efficient way of detecting stall

%Init
coaxial.state.forward_vel = 0;


iter_collective_axial = 10:10:60; %tweak this selection

axial_vel_range = 0:5:120;

CT_arr_axial = zeros(1,length(iter_collective_axial)*length(axial_vel_range));
CP_arr_axial = zeros(1,length(iter_collective_axial)*length(axial_vel_range));
axial_vel_arr = zeros(1,length(iter_collective_axial)*length(axial_vel_range));

method = 'leishman';

i = 1;
for collective_idx = 1:length(iter_collective_axial)
    coaxial.state.collective = iter_collective_axial(collective_idx);
    
    for axial_vel_idx = 1:length(axial_vel_range)
        
        coaxial.state.axial_vel = axial_vel_range(axial_vel_idx); %m/s
         
        %[coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        %coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT_axial(coaxial,atm,epsilon,plots,verbose);
        
        try %can be dangerous because it does not show obvious error messages
            tic
            [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,iter_collective_axial(collective_idx),'pitch_upper',method);
                                                                                    %trim(coaxial,atm,epsilon,CT_or_pitch,  trimvar,method)
            toc
        catch
            disp("Negative thrust")
            coaxial.state.CT = -1; %random negative number
        end
        if coaxial.state.CT < 0
            CT_arr_axial(i) = [];
            CP_arr_axial(i) = [];
            axial_vel_arr(i) = [];
            i = i-1;
        else
            CT_arr_axial(i) = coaxial.state.CT;
            CP_arr_axial(i) = coaxial.state.CP;
            axial_vel_arr(i) = axial_vel_range(axial_vel_idx);
        end
        i = i+1;
    end
end

figure(1)
scatter(axial_vel_arr,CT_arr_axial,'o')
figure(2)
scatter(axial_vel_arr,CP_arr_axial,'o')
%% Forward flight performance validation

%Init
warning('off','all')
coaxial.state.axial_vel = 0;
coaxial.state.side_vel = 0;
CT_desired = 0.0048/2; %THIS DEPENDS ON WHETHER IT IS COAXIAL OR SINGLE (COAXIAL = 2*SINGLE)
v_tip = 142.951; %m/s
coaxial.rotor(1).omega = v_tip/coaxial.rotor(1).R;

plots = false;
verbose = false;
debug = false;
method = 'leishman';

%Load data if current rotor is Harrington1

if ~strcmpi(coaxial.name,'Harrington1')
    data_mu = [];
    data_CP_brown_coax = [];
    data_CP_dingle_single = [];
else
    % Brown data coaxial
    data_CP_mu = readmatrix('H1_CP_mu_coax.csv'); %CP_mu
    data_mu = data_CP_mu(:,1);
    data_CP_brown_coax = data_CP_mu(:,2);
    
    %atm.rho = 1.3; %seemed to match better the plots from Brown and Dingledein 
    
    %Dingledein data coaxial
    data_HP_mu_coax_dingle = readmatrix('H1_CP_mu_coax_dingeldein.csv');
    data_CP_dingle_coax = data_HP_mu_coax_dingle(:,2)*745.7/(atm.rho*pi*coaxial.rotor(1).R^2*v_tip^3); %first conversion to W and then to dimensionless CP
    
    %Dingledein data single
    data_HP_mu_single_dingle = readmatrix('H1_CP_mu_single_dingeldein.csv');
    data_CP_dingle_single = data_HP_mu_single_dingle(:,2)*745.7/(atm.rho*pi*coaxial.rotor(1).R^2*v_tip^3); %first conversion to W and then to dimensionless CP
    
    %% Comparison of plots to check that the conversion is correct (units)
    %It is correct within a margin of error (which could come from WebplotDigitizer, etc)
    
    
%     hold on
%     scatter(data_mu,data_CP_brown_coax)
%     scatter(data_mu,data_CP_dingle_coax)
%     scatter(data_mu,data_CP_dingle_single)
    
end

advance_ratio_arr = linspace(0,0.35,15);
CP_arr = zeros(1,length(advance_ratio_arr));

i = 0;
for advance_ratio = advance_ratio_arr
    i=i+1;
    
    coaxial.state.forward_vel = advance_ratio*v_tip;
    
    [collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT",method);
   
    
    disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
    %disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
    %disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])
    

    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
    
    CP_arr(i) = coaxial.state.CP;   
    
end

% plots
hold on
scatter(data_mu,data_CP_brown_coax)
scatter(data_mu,data_CP_dingle_single)
plot(advance_ratio_arr,CP_arr)
xlabel('Advance ratio $\mu = \frac{V_T}{\Omega R}$','Interpreter','latex')
ylabel('$C_P$','interpreter','latex')
legend('Experiment Coaxial','Experiment Single','BEMT') % instead do set(gca,legend...) in the if statement above

warning('on','all')
%% Optimal collective plot single rotor at different axial speeds
%idea: instead of making this plot, just find the maximu CT/CP point and
%what collective it corresponds to. Then plot optimum collective vs axial
%flight speed and hopefully it is a horizontal line
axial_vel_range = 0:1:60;
iter_pitchdeg = 1:1:60;

coaxial.state.forward_vel = 0;

plots= false;
verbose= false;
debug = false;
method='leishman';

opt_collectives = zeros(1,length(axial_vel_range));

for vel_idx = 1:length(axial_vel_range)

    CT_arr = zeros(1,length(iter_pitchdeg));
    CP_arr = zeros(1,length(iter_pitchdeg));
    coaxial.state.axial_vel = axial_vel_range(vel_idx);
    for idx = 1:length(iter_pitchdeg)
            coaxial.state.collective = iter_pitchdeg(idx);
            try
                tic
                [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
                toc
                if CT < 0
                    CT_arr(idx) = nan;
                    CP_arr(idx) = nan;
                else
                    CT_arr(idx) = CT;
                    CP_arr(idx) = CP;
                end
            catch
                CT_arr(idx) = nan;
                CP_arr(idx) = nan;
                
            end
 
    end
    eff = CT_arr./CP_arr;
    [max_eff, collective_idx] = max(eff);
    opt_collectives(vel_idx) = iter_pitchdeg(collective_idx);
    %plot(CP_arr,CT_arr)
    %hold on
end

plot(axial_vel_range,opt_collectives)
xlabel('Axial velocity [m/s]')
ylabel('Optimum collective angle [deg]')
%% CT CP for different axial velocities

axial_vel_range = 0:1:60;
coaxial.state.forward_vel = 0;
coaxial.state.collective = 25; %collective in deg
method = 'leishman'; %airfoil took about 6 mins to run

for axial_vel_idx = 1:length(axial_vel_range)
    coaxial.state.axial_vel = axial_vel_range(axial_vel_idx); %m/s - comparison plots are for hover
    
    %%
    
    CT_arr = zeros(1,length(axial_vel_range));
    CP_arr = zeros(1,length(axial_vel_range));
    
    
    tic
    [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,coaxial.state.collective,'pitch_upper',method);
    toc
    if coaxial.state.CT < 0
        CT_arr(axial_vel_idx) = [];
        CP_arr(axial_vel_idx) = [];
        
    else
        CT_arr(axial_vel_idx) = coaxial.state.CT;
        CP_arr(axial_vel_idx) = coaxial.state.CP;

    end

    
    
end

figure(1)
%legend('-DynamicLegend');
hold all
set(gca,'FontSize',16)
plot(CP_arr,CT_arr,'LineWidth',2)

xlabel('$C_P$','Interpreter','latex')
ylabel('$C_T$','Interpreter','latex')

title(string(coaxial.name))
