%% Init

clear
close all
clc

addpath(genpath('./'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change parameters

coaxial.state.axial_vel = 2; %m/s 
coaxial.state.tangent_vel = 0; %m/s 
coaxial.state.trim= 1;
coaxial.state.collective = 20; %collective in deg

epsilon = 0.0001; %convergence accuracy for Fcf and lambda -> 0.0001

%warning('off')

%% Testing 

plots= true;
verbose= true;
debug = false;
method='leishman'; %'leishman','airfoil'

tic
[Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
toc


%%%%% GENERATE DIFFERENCE (between methods) DISK PLOTS SOMEWHERE!!!!!
%% Iteration to trim the coaxial rotor and produce CT-CP validation plots

iter_pitchdeg = 0:1:18;
method = 'leishman'; %airfoil took about 6 mins to run

coaxial.state.tangent_vel = 0; %m/s  - comparison plots are for hover
SMT = true;
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
coaxial.state.tangent_vel = 0; %m/s - hover
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

%SOMETHING SUPER WEIRD IS GOING ON HERE


iter_pitchdeg_axial = [20]; %tweak this selection

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

%% Forward flight performance validation

coaxial.state.axial_vel = 0;
CT_desired = 0.0048;

plots = false;
verbose = false;
debug = false;
method = 'leishman';

data_CP_mu = readmatrix('H1_CP_mu_coax.csv'); %

advance_ratio_arr = linspace(0,max(data_CP_mu(:,1)),15);
CP_arr = zeros(1,length(advance_ratio_arr));

i = 0;
for advance_ratio = advance_ratio_arr
    i=i+1;
    
    coaxial.state.tangent_vel = advance_ratio*coaxial.rotor(1).omega*coaxial.rotor(1).R;
    
    [collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"CT",method);
   
    
    disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(CT)])
    %disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
    %disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])
    
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);

    CP_arr(i) = coaxial.state.CP;   
    
end

%%
hold on
scatter(data_CP_mu(:,1),data_CP_mu(:,2),'o')
plot(advance_ratio_arr,CP_arr)
xlabel('Advance ratio $\mu = \frac{V_T}{\Omega R}$','Interpreter','latex')
ylabel('$C_P$','interpreter','latex')


%% Optimal collective plot single rotor at different axial speeds
%idea: instead of making this plot, just find the maximu CT/CP point and
%what collective it corresponds to. Then plot optimum collective vs axial
%flight speed and hopefully it is a horizontal line
axial_vel_range = 0:1:60;
iter_pitchdeg = 1:1:60;

coaxial.state.tangent_vel = 0;

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
coaxial.state.tangent_vel = 0;
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
