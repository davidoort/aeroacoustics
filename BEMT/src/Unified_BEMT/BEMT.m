function [Power, Forces, Moments, CT, CP, net_torque_coeff] = BEMT(rotorsystem,atm,epsilon,plots,verbose,method,debug)
%{

Power is a 1x2 matrix -> [P_upper; P_lower] [W]

Forces is a 3x2 matrix -> [Fx_upper, Fy_upper, Fz_upper; Fx_lower, Fy_lower, Fz_lower] [N]
using the right-handed coordinate system with x pointing forward, z pointing up

Moments is a 3x2 matrix -> [Mx_upper, My_upper, Mz_upper; Mx_lower, My_lower, Mz_lower] [Nm] 
using the same coordinate system as Forces

CT is the thrust coefficient 1x2 matrix -> [CT_u; CT_l];
CP is the torque coefficient 1x2 matrix -> [CP_u;CP_l];

net_torque_coeff = (CP_u-CP_l)/(CP_u+CP_l) [-] 

%}

%init is done inside the spinRotor function, but you have to specify which
%rotor (geometry) has to be spun and in which direction it has to be spun, with what
%collective (and cyclic later), lambda_P, lambda_T, method (leishman or airfoil)


%% Init 
if debug
    res_r = 10;
    res_psi = 20;
else
    res_r = 100;
    res_psi = 130; %careful of putting them equal
end

axial_vel = rotorsystem.state.axial_vel; %m/s
forward_vel = rotorsystem.state.forward_vel; %m/s
side_vel = rotorsystem.state.side_vel; %m/s

tangent_vel = norm([forward_vel,side_vel]); %m/s
sideslip = rotorsystem.state.sideslip(); %rad


%% Top Rotor

flowfield(1).lambda_P = axial_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %normalizing free stream axial velocity by tip velocity
flowfield(1).lambda_T = tangent_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %%normalizing free stream tangential velocity by tip velocity - FOR LATER

collective_u = rotorsystem.state.collective_u;
cyclic_u = [rotorsystem.state.cyclic_s, rotorsystem.state.cyclic_c];

spin_dir_u = 'CCW'; %should probably become a parameter in the rotor object
[Thrust_u, Torque_u, Power_u, CT_u, CP_u, dCT_u, dCP_u, lambda_u, Re_u, alpha_u,alpha_negatives_u, phi_u, ...
    F_u, weighted_swirl_ratio_u, FOM_u, velocity_dimensional_u, pitchdeg_u,r,dr,psi_u,chord_u] = spinRotor(rotorsystem.rotor(1),atm,spin_dir_u,collective_u,cyclic_u,flowfield(1).lambda_P,flowfield(1).lambda_T,method,epsilon);

net_torque_coeff = 0;

net_torque_dimensional = 0;

lambda_i_u_nans = lambda_u-flowfield(1).lambda_P;

lambda_i_u = lambda_i_u_nans;
lambda_i_u(isnan(lambda_i_u))=0;

%rotorsystem.state.CT = CT;
%rotorsystem.state.CP = CP;


% Calculation of forces and moments UPPER rotor

rho = atm.rho;
R_upper = rotorsystem.rotor(1).R;
omega_upper = rotorsystem.rotor(1).omega;

rotation_matrix = [cos(sideslip) sin(sideslip); -sin(sideslip) cos(sideslip)]; %derived visually in David's notebook

%create a general psi (which is positive in the CCW direction) - this
%fixes the problem of a non-sensical plot when upper rotor is CW
psi = abs(psi_u);


%Forces
%flow reference frame (x-axis pointing in the same direction as the free-stream velocity)
F_u_flow = [-sum(sum(dCP_u./r.*sin(psi))) * rho*pi*R_upper^2*omega_upper^2*R_upper^3/R_upper;... %x
             sum(sum(dCP_u./r.*cos(psi))) * rho*pi*R_upper^2*omega_upper^2*R_upper^3/R_upper]; %y
F_u_body = rotation_matrix*F_u_flow; %body reference frame

Fx_u = F_u_body(1);
Fy_u = F_u_body(2);
Fz_u = Thrust_u;

%Moments

%flow reference frame (x-axis pointing in the same direction as the free-stream velocity)
M_u_flow = [-sum(sum(dCT_u.*r.*sin(psi))) * rho*pi*R_upper^2*omega_upper^2*R_upper^2 * R_upper;... %x
             sum(sum(dCT_u.*r.*cos(psi))) * rho*pi*R_upper^2*omega_upper^2*R_upper^2 * R_upper]; %y - this should amount to zero for the upper rotor
M_u_body = rotation_matrix*M_u_flow; %body reference frame


Mx_u = M_u_body(1);
My_u = M_u_body(2);
Mz_u = (2*strcmpi(spin_dir_u,'CW')-1)*Torque_u;

%% Preliminary Bottom Rotor Output (there might not be a bottom rotor)

%Power
Power_l = 0;

%Forces
Fx_l = 0;
Fy_l = 0;
Fz_l = 0;
Thrust_l = 0; %the same as Fz_l

%Moments
Mx_l = 0;
My_l = 0;
Mz_l = 0; 
Torque_l = 0; %the same as Mz_l
%CT
CT_l = 0;

%CP
CP_l = 0;

%% Bottom Rotor
if strcmpi(rotorsystem.type,"coaxial")
    %% Init
    
    collective_l = rotorsystem.state.collective_l;
    cyclic_l = -[rotorsystem.state.cyclic_s, rotorsystem.state.cyclic_c];
    
    flowfield(2).lambda_P = axial_vel/(rotorsystem.rotor(2).omega*rotorsystem.rotor(2).R)*ones(res_psi,res_r); %normalizing free stream axial velocity by tip velocity
    flowfield(2).lambda_T = tangent_vel/(rotorsystem.rotor(2).omega*rotorsystem.rotor(2).R)*ones(res_psi,res_r); %%normalizing free stream tangential velocity by tip velocity - FOR LATER
   
    
    %% Calculate inflow for bottom rotor
    
    separation = rotorsystem.params.interrotor_spacing*rotorsystem.rotor(1).R;
    
    %radial_contraction = getContraction(separation,CT) %Landgrebe
    
    radial_contraction = rotorsystem.params.rd;
    
    mean_downwash = mean(mean(lambda_i_u))*rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R; 
    %check if in hover to avoid singularities when lambda_i =0
    if axial_vel==0 && forward_vel ==0
        skew_angle = 0; %for any lambda_i in hover
    else
        skew_angle = atan(flowfield(1).lambda_T./(flowfield(1).lambda_P+lambda_i_u)); %more accurate (disk element by disk element) - gives non circular marks
    end
    r_contracted = radial_contraction*r;
    psi_contracted= psi_u;
    lambda_contracted = lambda_i_u/(radial_contraction)^2;
    
    x_contracted = r_contracted.*cos(psi_contracted);
    y_contracted = r_contracted.*sin(psi_contracted);
    
    x_skewed = tan(skew_angle)*separation+x_contracted;
    y_skewed = y_contracted;
    
    %lambda_contracted(r_contracted<radial_contraction*rotorsystem.rotor(1).hub_radial_fraction) = 0;
    
    lambda_skewed = lambda_contracted; %the basis changes
    
    
    
    r_skewed = sqrt(x_skewed.^2+y_skewed.^2);
    a = atan2(y_skewed,x_skewed); 
    psi_skewed = a .* (a >= 0) + (a + 2 * pi) .* (a < 0); %now psi still goes from 0 to 2pi (will be handy for interpolation)
    
   
    % Trying to solve the problem in cartesian
    %{
    [X,Y] = meshgrid(-1:0.01:1);
    
    lambda_i_nan = griddata(x_skewed,y_skewed,lambda_skewed,X,Y);
    
    lambda_i = lambda_i_nan;
    lambda_i(isnan(lambda_i_nan))=0;
    
    plot3(x_skewed,y_skewed,lambda_skewed,'o')
    hold on
    %mesh(X,Y,lambda_i)
    
    R = sqrt(X.^2+Y.^2);
    A = atan2(Y,X); 
    PSI = A .* (A >= 0) + (A + 2 * pi) .* (A < 0); %now psi still goes from 0 to 2pi (will be handy for interpolation)

    %mesh(R.*cos(PSI),R.*sin(PSI),lambda_i)
    
    X(R>1) = nan;
    X(R<rotorsystem.rotor(1).hub_radial_fraction) = nan;
    Y(R>1) = nan;
    Y(R<rotorsystem.rotor(1).hub_radial_fraction) = nan;
    lambda_i(R>1) = nan;
    lambda_i(R<rotorsystem.rotor(1).hub_radial_fraction) = nan;
    
    R = sqrt(X.^2+Y.^2);
    A = atan2(Y,X); 
    PSI = A .* (A >= 0) + (A + 2 * pi) .* (A < 0); %now psi still goes from 0 to 2pi (will be handy for interpolation)

    mesh(R.*cos(PSI),R.*sin(PSI),lambda_i)
    %}
    
    
    not_outside_of_disk = r_skewed<1;
    not_inside_of_hub = r_skewed>rotorsystem.rotor(1).hub_radial_fraction;
    inside_disk = not_outside_of_disk.*not_inside_of_hub;
    lambda_skewed= lambda_skewed.*(inside_disk);
    

    %Other interpolation methods
    %{
    %lambda_P_bottom = interp2(r_skewed,psi_skewed,lambda_skewed,r,psi,'linear',flowfield(2).lambda_P(1,1)); 
%     interpolant = scatteredInterpolant(r_skewed', psi_skewed', lambda_skewed');
%     lambda_P_bottom_nans = interpolant(r, psi);
%     lambda_P_bottom = lambda_P_bottom_nans;
%     lambda_P_bottom(isnan(lambda_P_bottom_nans)) = 0; 
    %} 
    
    % Interpolation in cartesian coordinates (to avoid extrapolation)
    
    lambda_P_bottom_ind_nans = griddata(r_skewed.*cos(psi_skewed),r_skewed.*sin(psi_skewed),lambda_skewed,r.*cos(psi_u),r.*sin(psi_u),'cubic'); 
    
    %lambda_P_bottom_ind_nans = griddata(r_skewed,psi_skewed,lambda_skewed,r,psi_u,'cubic'); 
    lambda_P_bottom_ind = lambda_P_bottom_ind_nans;
    lambda_P_bottom_ind(isnan(lambda_P_bottom_ind_nans)) = 0; 
    
    %Debugging
    %{
    figure(1)
    plot3(r_skewed.*cos(psi_skewed),r_skewed.*sin(psi_skewed),lambda_skewed,'o')
    hold on
    mesh(r.*cos(psi_u),r.*sin(psi_u),lambda_P_bottom_ind)
    
    figure(2) %let's check if it is actually extrapolating
    plot3(r_skewed,psi_skewed,lambda_skewed,'o')
    hold on
    mesh(r,psi_u,lambda_P_bottom_ind)
    %}
    
    % Add the free-stream axial velocity of the rotorsystem
    lambda_P_bottom = lambda_P_bottom_ind + flowfield(2).lambda_P(1,1)*ones(size(lambda_P_bottom_ind)); %lambda_P (freestream) as extrapval
    
    
    %% Bottom Rotor
    
    spin_dir_l = 'CW';
    [Thrust_l, Torque_l, Power_l, CT_l, CP_l, dCT_l, dCP_l, lambda_l, Re_l, alpha_l,alpha_negatives_l, phi_l, ...
    F_l, weighted_swirl_ratio_l, FOM_l, velocity_dimensional_l,pitchdeg_l,r,dr,psi_l,chord_l] = spinRotor(rotorsystem.rotor(2),atm,spin_dir_l,collective_l,cyclic_l,lambda_P_bottom,flowfield(2).lambda_T,method,epsilon);

    %% Output Vars

        
    %Thrust = Thrust_u+Thrust_l;
    %Torque = Torque_u+Torque_l;
    %Power = Power_u+Power_l;
    %CT = CT_u+CT_l;
    %CP = CP_u+CP_l;
    
    
    % Calculation of forces and moments UPPER rotor
    
    R_lower = rotorsystem.rotor(2).R;
    omega_lower = rotorsystem.rotor(2).omega;
    
    %Forces
    %flow reference frame (x-axis pointing in the same direction as the free-stream velocity)
    F_l_flow = [-sum(sum(dCP_l./r.*sin(psi))) * rho*pi*R_lower^2*omega_lower^2*R_lower^3/R_lower;... %x
        sum(sum(dCP_l./r.*cos(psi))) * rho*pi*R_lower^2*omega_lower^2*R_lower^3/R_lower]; %y
    F_l_body = rotation_matrix*F_l_flow; %body reference frame
    
    Fx_l = F_l_body(1);
    Fy_l = F_l_body(2);
    Fz_l = Thrust_l;
    
    %Moments
    
    %flow reference frame (x-axis pointing in the same direction as the free-stream velocity)
    M_l_flow = [-sum(sum(dCT_l.*r.*sin(psi))) * rho*pi*R_lower^2*omega_lower^2*R_lower^2 * R_lower;... %x
        sum(sum(dCT_l.*r.*cos(psi))) * rho*pi*R_lower^2*omega_lower^2*R_lower^2 * R_lower]; %y - this should amount to zero for the upper rotor
    M_l_body = rotation_matrix*M_l_flow; %body reference frame
    
    
    Mx_l = M_l_body(1);
    My_l = M_l_body(2);
    Mz_l = (2*strcmpi(spin_dir_l,'CW')-1)*Torque_l;
    
    
    FOM_coax = rotorsystem.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper
    
    net_torque_coeff = (CP_u-CP_l)/(CP_u+CP_l);
    net_torque_dimensional = Torque_u-Torque_l;
    
    %% Noise - strictly for coaxial rotors
    
    %since I am doing the implementation for a coaxial rotor system
    
    %this will return the loading harmonics up to a high number k which can
    %then be truncated in hanson_aeroacoustics
    
    [CLk_u,CDk_u] = getLoadingHarmonics(rotorsystem.rotor,dCT_u,dCP_u,phi_u,velocity_dimensional_u,chord_u,r,psi_u);
    [CLk_l,CDk_l] = getLoadingHarmonics(rotorsystem.rotor,dCT_l,dCP_l,phi_l,velocity_dimensional_l,chord_l,r,psi_l);
    
    CLk.upper = CLk_u;
    CLk.lower = CLk_l;
    CDk.upper = CDk_u;
    CDk.lower = CDk_l;
    
    %hanson_acoustics is simply put the code version of equation 1 in
    %Hanson's paper
    [dB] = hanson_acoustics(coaxial,atm,r,dr,CLk,CDk);

    
end

     



%% Outputs

%Power
Power = [Power_u;Power_l];

%Forces
Forces = [Fx_u, Fy_u, Fz_u; Fx_l, Fy_l, Fz_l];

%Moments
Moments = [Mx_u, My_u, Mz_u; Mx_l, My_l, Mz_l];

%CT
CT = [CT_u; CT_l];

%CP
CP = [CP_u;CP_l];




%% Display output text

if verbose
    disp('                                                  ')
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' rotor specs:'))
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' power coefficient [-] ',num2str(CP_u+CP_l)))
    disp(strcat(rotorsystem.type,' thrust coefficient [-] ',num2str(CT_u+CT_l)))
    disp(['Total thrust [N]',' ',num2str(Thrust_u+Thrust_l)])
    disp(['Total power [W]',' ',num2str(Power_u+Power_l)])
    disp(['Total torque [Nm]',' ',num2str(Torque_u+Torque_l)])
    disp(['Omega [rad/s]',' ',num2str(rotorsystem.rotor(1).omega)])
    disp(['Swirl upper rotor [rad/s] ',num2str(weighted_swirl_ratio_u*rotorsystem.rotor(1).omega)])
    disp(['Max/Min AoA upper rotor [deg] ',num2str(max(max(alpha_negatives_u))),' / ',num2str(min(min(alpha_negatives_u)))])
    disp(['Max/Min/0.7R Re number upper [-] ',num2str(max(max(Re_u))),' / ',num2str(min(min(Re_u))),' / ',num2str(interp1(r(1,:),Re_u(1,:),0.7))])
    disp(['Pitch upper rotor [deg] ',num2str(pitchdeg_u(1,:))])        
    if strcmpi(rotorsystem.type,"coaxial")
        disp(['Net torque coefficient (u-l)/l [-]',' ',num2str(net_torque_coeff)])
        disp(['Net torque (u-l) [Nm]',' ',num2str(net_torque_dimensional)])
        disp(['Collective upper/lower rotor [deg] ',num2str(collective_u),' / ',num2str(collective_l)])
        disp(['Max/Min AoA lower rotor [deg] ',num2str(max(max(alpha_negatives_l))),' / ',num2str(min(min(alpha_negatives_l)))])
        disp(['Max/Min/0.7R Re number lower [-] ',num2str(max(max(Re_l))),' / ',num2str(min(min(Re_l))),' / ',num2str(interp1(r(1,:),Re_l(1,:),0.7))])
        disp(strcat(rotorsystem.type,' FOM [-] ',num2str(FOM_coax)))
    else
        disp(strcat(rotorsystem.type,' FOM [-] ',num2str(FOM_u)))
    end
    
end   


%% Plots

if plots
    %% DISK PLOTS
    figure(99); clf;
    title('Upper Rotor')
    subplot(2, 3, 1)
    hold on
    diskPlot(r,psi,alpha_u,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
    subplot(2,3,2)
    diskPlot(r,psi,dCT_u,rotorsystem.state,'$dC_T$ [-]')
    subplot(2,3,3)
    diskPlot(r,psi,dCP_u,rotorsystem.state,'$dC_P$ [-]')
    subplot(2,3,4)
    diskPlot(r,psi,lambda_u,rotorsystem.state, '$\lambda_u = \frac{v_u+V_P}{v_{tip}}$ [-]')
    subplot(2,3,5)
    diskPlot(r,psi,F_u,rotorsystem.state,'Prandtl [-]')
    subplot(2,3,6)
    diskPlot(r,psi,rad2deg(phi_u),rotorsystem.state,'$\phi$ [deg]')
    
    %For leonardo times
%     figure(101); clf;
%     title('Upper Rotor')
%     subplot(2, 2, 1)
%     hold on
%     diskPlot(r,psi,alpha_u,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
%     subplot(2,2,2)
%     diskPlot(r,psi,dCT_u,rotorsystem.state,'Thrust Coefficient [-]')
%     subplot(2,2,3)
%     diskPlot(r,psi,dCP_u,rotorsystem.state,'Power Coefficient [-]')
%     subplot(2,2,4)
%     diskPlot(r,psi,lambda_u,rotorsystem.state, 'Inflow ratio $\frac{v_u+V_P}{v_{tip}}$ [-]')

    if strcmpi(rotorsystem.type,"coaxial")
        figure(100); clf;
        title('Lower Rotor')
        subplot(2, 3, 1)
        hold on
        diskPlot(r,psi,alpha_l,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
        subplot(2,3,2)
        diskPlot(r,psi,dCT_l,rotorsystem.state,'$dC_T$ [-]')
        subplot(2,3,3)
        diskPlot(r,psi,dCP_l,rotorsystem.state,'$dC_P$ [-]')
        subplot(2,3,4)
        diskPlot(r,psi,lambda_l,rotorsystem.state,'$\lambda_l = \frac{v_l+V_P+v_u/a^2}{v_{tip}}$ [-]')
        subplot(2,3,5)
        diskPlot(r,psi,F_l,rotorsystem.state,'Prandtl [-]')
        subplot(2,3,6)
        diskPlot(r,psi,rad2deg(phi_l),rotorsystem.state,'$\phi$ [-]')

        
        %for leonardo times
%         figure(102); clf;
%         title('Lower Rotor')
%         subplot(2, 2, 1)
%         hold on
%         diskPlot(r,psi,alpha_l,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
%         subplot(2,2,2)
%         diskPlot(r,psi,dCT_l,rotorsystem.state,'Thrust Coefficient [-]')
%         subplot(2,2,3)
%         diskPlot(r,psi,dCP_l,rotorsystem.state,'Power Coefficient [-]')
%         subplot(2,2,4)
%         diskPlot(r,psi,lambda_l,rotorsystem.state,'Inflow ratio $\frac{v_l+V_P+v_u/a^2}{v_{tip}}$ [-]')
         
    end    
    %% Validation plots
    if axial_vel==0 && tangent_vel==0
        if abs(sum(CT)-0.004)<0.00001 && net_torque_dimensional<0.1 && strcmpi(rotorsystem.name,"Harrington1")  
            data_inflow_upper = readmatrix('H1_inflow_FVM_fig10a.csv');
            r1 = data_inflow_upper(:,1); lambda1 = data_inflow_upper(:,2);
            
            data_dCT_upper = readmatrix('H1_dCT_FVM_fig11a.csv');
            r3 = data_dCT_upper(:,1); dCT1 = data_dCT_upper(:,2);
            
            data_dCP_upper = readmatrix('H1_dCP_FVM_fig12a.csv');
            r5 = data_dCP_upper(:,1); dCP1 = data_dCP_upper(:,2);
            
            data_inflow_lower = readmatrix('H1_inflow_FVM_fig10b.csv');
            r2 = data_inflow_lower(:,1); lambda2 = data_inflow_lower(:,2);
            
            data_dCT_lower = readmatrix('H1_dCT_FVM_fig11b.csv');
            r4 = data_dCT_lower(:,1); dCT2 = data_dCT_lower(:,2);
            
            data_dCP_lower = readmatrix('H1_dCP_FVM_fig12b.csv');
            r6 = data_dCP_lower(:,1); dCP2 = data_dCP_lower(:,2);

        else
            
            r1=[];r2=[];r3=[];r4=[];r5=[];r6=[];
            lambda1=[];lambda2=[];dCT1=[];dCT2=[];dCP1=[];dCP2=[];
            
        end
        
        %%%%%
        if strcmpi(rotorsystem.type,"coaxial")
            figure(1); clf;
            subplot(2, 3, 1)
            hold on
            scatter(r1,lambda1)
            plot(r, lambda_u(1,:), 'b-.')
            title('Inflow ratio vs radius - Top')
            xlabel('r/R')
            ylabel('Inflow ratio')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(2, 3, 2)
            plot(r, F_u(1,:), 'b-.')
            title('Prandtl tip loss vs radius - Top')
            xlabel('r/R')
            ylabel('Prandtl tip loss')

            subplot(2, 3, 3)
            plot(r, alpha_u(1,:), 'b-.')
            title('AoA vs radius - Top')
            xlabel('r/R')
            ylabel('$\alpha$ [deg]','interpreter','latex')
            ylim([-20 Inf])

            subplot(2, 3, 4)
            hold on
            scatter(r2,lambda2)
            plot(r, lambda_l(1,:), 'b-.')
            title('Inflow ratio vs radius - Bottom')
            xlabel('r/R')
            ylabel('Inflow ratio')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(2, 3, 5)
            plot(r, F_l(1,:), 'b-.')
            title('Prandtl tip loss vs radius - Bottom')
            xlabel('r/R')
            ylabel('Prandtl tip loss')

            subplot(2, 3, 6)
            plot(r, alpha_l(1,:), 'b-.')
            title('AoA vs radius - Bottom')
            xlabel('r/R')
            ylabel('$\alpha$ [deg]','interpreter','latex')
            ylim([-20 Inf])

            %%%%%%%

            figure(2); clf;
            subplot(2, 2, 1)
            hold on
            scatter(r3,dCT1)
            plot(r, sum(dCT_u)/dr, 'b-.')
            title('dCTu vs radius - Top')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCTu/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(2, 2, 2)
            hold on
            scatter(r5,dCP1)
            plot(r, sum(dCP_u)/dr, 'b-.')
            title('dCpu vs radius - Top')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCpu/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(2, 2, 3)
            hold on
            scatter(r4,dCT2)
            plot(r, sum(dCT_l)/dr, 'b-.')
            title('dCTl vs radius - Bottom')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCtl/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(2, 2, 4)
            hold on
            scatter(r6,dCP2)
            plot(r, sum(dCP_l)/dr, 'b-.')
            title('dCpl vs radius - Bottom')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCpl/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')
        
        elseif strcmpi(rotorsystem.type,"single")
            figure(1); clf;
            subplot(1, 3, 1)
            hold on
            scatter(r1,lambda1)
            plot(r, lambda_u(1,:), 'b-.')
            title('Inflow ratio vs radius - Top')
            xlabel('r/R')
            ylabel('Inflow ratio')
            xlim([0.2 1])
            legend('FVM','BEMT')

           
            subplot(1, 3, 2)
            hold on
            scatter(r3,dCT1)
            plot(r, sum(dCT_u)/dr, 'b-.')
            title('dCTu vs radius - Top')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCTu/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')

            subplot(1, 3, 3)
            hold on
            scatter(r5,dCP1)
            plot(r, sum(dCP_u)/dr, 'b-.')
            title('dCpu vs radius - Top')
            xlabel('r/R')
            ylabel('Non-dimensionalized dCpu/dr')
            xlim([0.2 1])
            legend('FVM','BEMT')
        end
    end
end

    
end






