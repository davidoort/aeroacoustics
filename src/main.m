%% Init

clear
close all
clc

addpath(genpath('../'))

% Instantiate objects

atm = Atmosphere();
coaxial = Rotor();

% Change flight parameters

coaxial.state.axial_vel = 0; %m/s 
coaxial.state.forward_vel = 0; %m/s 
coaxial.state.side_vel = 0; %m/s

% Control Inputs
%coaxial.rotor(1).omega = 125;
coaxial.state.collective_u = 60; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.collective_l = 60; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.cyclic_s = 0; %sine term for cyclic (gets multiplied by sin(azimuth))
coaxial.state.cyclic_c = 0; %cosine term for cyclic (gets multiplied by cos(azimuth))

epsilon = 0.0001; %convergence accuracy for Prandtl tip function and inflow ratio
timelimit = inf;
acoustics = false;


%% Testing 

plots= true;
verbose= true;
debug = false;
method='leishman'; %'leishman' (fast), 'airfoil' (slow and less robust convergence but potentially more accurate)
acoustics = true;

[Power, Forces, Moments, CT, CP, net_torque_coeff,sound] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);
 

%% Axial flight plots

close all

warning('off','all')

coaxial.state.forward_vel = 0;
coaxial.state.side_vel = 0;

method = 'leishman';
plots = false;
verbose = false;
debug = false;
%timelimit = 4; %seconds - as measured by tic toc

calculations_per_pitch = 20;

%SET UP COAXIAL AND SINGLE ROTOR PARAMETERS AND LOAD VALIDATION DATA
if strcmpi(coaxial.name,'NACA_single')
    if strcmpi(coaxial.rotor(1).twist_type,'linear')
        twist_increase = coaxial.rotor(1).twistdeg;
    elseif strcmpi(coaxial.rotor(1).twist_type,'ideal')
        twist_increase = coaxial.rotor(1).pitch_root;
    end
    
    iter_collective_axial = twist_increase*ones(1,7)+[22 26 30 35 40 45 50]; %tweak this selection - 7 lines should match pretty well
    
    %For more compact and higher resolution plots (and to solve the problem of blank regions in CP)
    %                   lb   ub
    axial_vel_range_arr = [25 75;
                           40 90;
                           60 110;
                           75 140;
                           100 160
                           125 190
                           160 240];
    
    %Load validation data
    v_tip = coaxial.rotor(1).R*coaxial.rotor(1).omega; %for de-normalization of axes
    solidity = coaxial.rotor(1).solidity;

    data_axial_flight_CP = readmatrix('NACA_single_CP_axial.csv');
    data_axial_flight_CT = readmatrix('NACA_single_CT_axial.csv');
    
    EXP_axial_vel_cp = data_axial_flight_CP(:,1)*v_tip; %m/s
    EXP_axial_vel_ct = data_axial_flight_CT(:,1)*v_tip; %m/s
    EXP_CT_axial = data_axial_flight_CT(:,2)*solidity;
    EXP_CP_axial = data_axial_flight_CP(:,2)*solidity;
    
elseif strcmpi(coaxial.name,'NACA_coax')
    
    %The following if statement is assuming that the rotors are not
    %different in terms of blades they use.
    
    if strcmpi(coaxial.rotor(1).twist_type,'linear')
        twist_increase = coaxial.rotor(1).twistdeg;
    elseif strcmpi(coaxial.rotor(1).twist_type,'ideal')
        twist_increase = coaxial.rotor(1).pitch_root;
    end
    
    
    %axial_vel_range_arr = repmat([0 105],[7,1]);
    
    axial_vel_range_arr =  [0 30;
                            0 35;
                            0 45;
                            10 55;
                            10 65
                            10 75
                            15 110
                            15 110];

    %UNKNOWN PARAMETER TWEAKING
    iter_collective_axial = twist_increase*ones(1,8)+[18 22 26 30 37 43 48 54]; %tweak this selection - 7 lines should match pretty well


    %Load validation data
    v_tip = coaxial.rotor(1).R*coaxial.rotor(1).omega; %for de-normalization of axes
    solidity = coaxial.rotor(1).solidity;
    data_axial_flight_CP = readmatrix('NACA_coaxial_CP_axial.csv');
    data_axial_flight_CT = readmatrix('NACA_coaxial_CT_axial.csv');
    
    EXP_axial_vel_cp = data_axial_flight_CP(:,1)*v_tip; %m/s
    EXP_axial_vel_ct = data_axial_flight_CT(:,1)*v_tip; %m/s
    EXP_CT_axial = data_axial_flight_CT(:,2)*solidity;
    EXP_CP_axial = data_axial_flight_CP(:,2)*solidity;
    
else
    error('Select a NACA propeller in Rotor')
    
end

%Iteration
for collective_idx = 1:length(iter_collective_axial)
    coaxial.state.collective = iter_collective_axial(collective_idx);
    i = 1;
    disp(['Spinning at collective setting of ', num2str(coaxial.state.collective), ' deg'])
    axial_vel_range = linspace(axial_vel_range_arr(collective_idx,1),axial_vel_range_arr(collective_idx,2),calculations_per_pitch);
    
    %axial_vel_range = 0:5:200;
    
    CT_arr_axial = zeros(1,length(axial_vel_range));
    CP_arr_axial = zeros(1,length(axial_vel_range));
    axial_vel_arr = zeros(1,length(axial_vel_range));
    
    for axial_vel_idx = 1:length(axial_vel_range)
        
        coaxial.state.axial_vel = axial_vel_range(axial_vel_idx); %m/s
         
        try %can be dangerous because it does not show obvious error messages
            
            [~, ~, ~, ~] = trim(coaxial,atm,epsilon,iter_collective_axial(collective_idx),'yaw',method,timelimit);
            [~, ~, ~, CT, CP, ~] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);
                
        catch
            disp("Negative thrust")
            CT = -1; %random negative number
        end
        
        toc
        
        if sum(CT) < 0 || (strcmpi(coaxial.name,'NACA_single') && sum(CT>0.07)) || toc>=timelimit 
            CT_arr_axial(i) = [];
            CP_arr_axial(i) = [];
            axial_vel_arr(i) = [];
            i = i-1;
            text_out = ' interrupted trim routine';
        else
            CT_arr_axial(i) = sum(CT);
            CP_arr_axial(i) = sum(CP);
            axial_vel_arr(i) = axial_vel_range(axial_vel_idx);
            text_out = '';
        end
        
        disp([toc, text_out])
        
        i = i+1;
    end
    
    figure(1)
    BEMT_ct = plot(axial_vel_arr,CT_arr_axial,'Color','b');
    hold on
    
    figure(2)
    BEMT_cp = plot(axial_vel_arr,CP_arr_axial,'Color','b');
    hold on
    
end

%PLOTS
figure(1)
hold on
experiment = scatter(EXP_axial_vel_ct, EXP_CT_axial, 'o');
legend([BEMT_ct(1), experiment(1)], 'BEMT', 'Experiment')
xlabel('Axial velocity $V_P$ [m/s]','Interpreter','Latex')
ylabel('$C_T$ [-]','Interpreter','Latex')

figure(2)
hold on
experiment = scatter(EXP_axial_vel_cp, EXP_CP_axial, 'o');
legend([BEMT_cp(1), experiment(1)], 'BEMT', 'Experiment')
xlabel('Axial velocity $V_P$ [m/s]','Interpreter','Latex')
ylabel('$C_P$ [-]','Interpreter','Latex')

warning('on','all')

%% Forward flight performance validation

%Init
warning('off','all')
coaxial.state.axial_vel = 0;
coaxial.state.side_vel = 0;

if strcmpi(coaxial.type,'single')
    CT_desired = 0.0048/2; 
elseif strcmpi(coaxial.type,'coaxial')
    CT_desired = 0.0048; 
end

v_tip = 142.951; %m/s
coaxial.rotor(1).omega = v_tip/coaxial.rotor(1).R;

plots = false;
verbose = false;
debug = false;
method = 'leishman';
acoustics = false;
%timelimit = 10;
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
    
end

advance_ratio_arr = linspace(0,0.35,15);
CP_arr = zeros(1,length(advance_ratio_arr));

i = 0;
for advance_ratio = advance_ratio_arr
    i=i+1;
    
    coaxial.state.forward_vel = advance_ratio*v_tip;
    
    tic
    [collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"thrust",method,timelimit);
    toc
    
    disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(sum(CT))])
    

    [Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);
    
    
    % For this small coaxial helicopter, an equivalent flat-plate parasite-drag area of 10 square feet was used. 
    A = 0.92903; %10sq->m^2
    Cd = 2; 
    
    %parasite_power = F*V = 1/2 * rho * V^3 * A *Cd
    parasite_power = 0.5*atm.rho*coaxial.state.forward_vel^3*A*Cd;
    
    CP_arr(i) = sum(CP) + parasite_power/(atm.rho*pi*coaxial.rotor(1).R^5*coaxial.rotor(1).omega^3);   
    
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

%% Iteration to trim the coaxial rotor and produce CT-CP validation plots

warning('off','all')

iter_pitchdeg = 0:1:18;
method = 'leishman'; %airfoil took about 6 mins to run
verbose = false;
plots = false;
debug = true;

coaxial.state.forward_vel = 0; %m/s  - comparison plots are for hover
coaxial.state.axial_vel = 0; %m/s - comparison plots are for hover
coaxial.state.side_vel = 0; %m/s - comparison plots are for hover
SMT = true;

rotors = ["single","coaxial"];

for rotor_type = rotors
    
    CT_arr = zeros(1,length(iter_pitchdeg));
    CP_arr = zeros(1,length(iter_pitchdeg));
    %C_T_SMT_arr = zeros(1,length(iter_pitchdeg));
    C_P_SMT_arr = zeros(1,length(iter_pitchdeg));
    
    coaxial.type = rotor_type;
    
    for idx = 1:length(iter_pitchdeg)
        tic
        [~, ~, ~, ~] = trim(coaxial,atm,epsilon,iter_pitchdeg(idx),'yaw',method,timelimit);
        toc
        if coaxial.state.CT < 0
            CT_arr(idx) = [];
            CP_arr(idx) = [];
            
        else
            CT_arr(idx) = sum(coaxial.state.CT);
            CP_arr(idx) = sum(coaxial.state.CP);
            
            if strcmpi(coaxial.type,"coaxial")
                % SMT Coaxial
                C_P_SMT_arr(idx) = 0.5*coaxial.params.kappaint*coaxial.params.kappa*(sum(coaxial.state.CT))^(3/2)+(coaxial.rotor(2).solidity + coaxial.rotor(1).solidity)*coaxial.rotor(1).aero.Cd0/8; %using Cd0 of bottom rotor but could be for top rotor, this is for validation purposes
            else
                %SMT Single
                C_P_SMT_arr(idx) = coaxial.params.kappa*sum(coaxial.state.CT)^(3/2)/sqrt(2) + coaxial.rotor(1).solidity*coaxial.rotor(1).aero.Cd0/8;
                
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
        CP_exp = val_data(:,1);
        CT_exp = val_data(:,2);
    elseif strcmpi(coaxial.name,"Harrington2")
        ax_xlim = [0 0.001];
        ax_ylim = [0 0.01];
        if strcmpi(coaxial.type,"single")
            val_data =readmatrix('H2_single_fig3.csv');
        else
            val_data =readmatrix('H2_coax_fig3.csv');
        end
        CP_exp = val_data(:,1);
        CT_exp = val_data(:,2);
    else
        ax_xlim = [0 Inf];
        ax_ylim = [0 Inf];
        CP_exp = [];
        CT_exp = [];
    end

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

warning('on','all')

%% Verification plots with FVM for inflow, CT and CP distributions. Trim the coaxial rotor at a specified thrust coefficient

CT_desired = 0.004;
method='leishman';

coaxial.state.collective_u = 20; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.collective_l = 20; %UPPER rotor collective in deg - geometric pitch angle at the root of the UPPER rotor blades!
coaxial.state.cyclic_s = 0; %sine term for cyclic (gets multiplied by sin(azimuth))
coaxial.state.cyclic_c = 0; %cosine term for cyclic (gets multiplied by cos(azimuth))

coaxial.state.axial_vel = 0; %m/s - hover
coaxial.state.forward_vel = 0; %m/s - hover
coaxial.state.side_vel = 0; %m/s - hover
epsilon = 0.001; %convergence accuracy for Fcf and lambda -> 0.0001

[collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_desired,"thrust",method,timelimit);

disp(['Converged to ',num2str(net_torque_dimensional),' net torque [Nm] and CT = ',num2str(sum(CT))])
disp(['Pitch upper rotor = ', num2str(collective_u),' deg'])
disp(['Pitch lower rotor = ', num2str(collective_l),' deg'])


plots = true;
verbose = false;
debug = false;

[Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);

%% Optimal collective plot single rotor at different axial speeds (slow)
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
                [Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);
                toc
                if CT < 0
                    CT_arr(idx) = nan;
                    CP_arr(idx) = nan;
                else
                    CT_arr(idx) = sum(CT);
                    CP_arr(idx) = sum(CP);
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

%% CT CP for different axial velocities (slow)

axial_vel_range = 0:1:60;
coaxial.state.forward_vel = 0;
coaxial.state.collective = 25; %collective in deg
method = 'leishman'; %airfoil took about 6 mins to run

for axial_vel_idx = 1:length(axial_vel_range)
    coaxial.state.axial_vel = axial_vel_range(axial_vel_idx); %m/s - comparison plots are for hover
    
    CT_arr = zeros(1,length(axial_vel_range));
    CP_arr = zeros(1,length(axial_vel_range));
    
    tic
    [collective_u, collective_l, net_torque_dimensional, coaxial.state.CT] = trim(coaxial,atm,epsilon,coaxial.state.collective,'yaw',method,timelimit);
    toc
    if coaxial.state.CT < 0
        CT_arr(axial_vel_idx) = [];
        CP_arr(axial_vel_idx) = [];
        
    else
        CT_arr(axial_vel_idx) = sum(coaxial.state.CT);
        CP_arr(axial_vel_idx) = sum(coaxial.state.CP);

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
