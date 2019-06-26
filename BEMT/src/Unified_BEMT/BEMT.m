function [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT(rotorsystem,atm,epsilon,plots,verbose,method,debug)


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

axial_vel = rotorsystem.state.axial_vel;
tangent_vel = rotorsystem.state.tangent_vel;

%% Sanity check - should allow ot do try/catch instead of going into infinite loops

if strcmpi(method,'leishman') && (rotorsystem.rotor(1).pitch_root>rotorsystem.state.collective && strcmpi(rotorsystem.rotor(1).twist_type,'ideal')) || (rotorsystem.rotor(1).twistdeg>rotorsystem.state.collective && strcmpi(rotorsystem.rotor(1).twist_type,'linear'))
    error('Rotor upper: Your twist is larger than your collective! Negative pitch angles near the tip will break Leishman')
end

%% Top Rotor

flowfield(1).lambda_P = axial_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %normalizing free stream axial velocity by tip velocity
flowfield(1).lambda_T = tangent_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %%normalizing free stream tangential velocity by tip velocity - FOR LATER

collective_u = rotorsystem.state.collective;

[Thrust, Torque, Power, CT, CP, dCT_u, dCP_u, lambda_u, Re_u, alpha_u,alpha_negatives_u, phi_u, ...
    F_u, weighted_swirl_ratio_u, FOM_u, velocity_dimensional_u,pitchdeg_u,r,dr,psi_u] = spinRotor(rotorsystem.rotor(1),atm,'CCW',collective_u,flowfield(1).lambda_P,flowfield(1).lambda_T,method,epsilon);

net_torque_coeff = 0;

net_torque_dimensional = 0;

lambda_i_u_nans = lambda_u-flowfield(1).lambda_P;

lambda_i_u = lambda_i_u_nans;
lambda_i_u(isnan(lambda_i_u))=0;

if strcmpi(rotorsystem.type,"coaxial")
    %% Init
    
    Thrust_u = Thrust;
    Torque_u = Torque;
    Power_u = Power;
    CT_u = CT;
    CP_u = CP;

    collective_l = collective_u*rotorsystem.state.trim;
    
    flowfield(2).lambda_P = axial_vel/(rotorsystem.rotor(2).omega*rotorsystem.rotor(2).R)*ones(res_psi,res_r); %normalizing free stream axial velocity by tip velocity
    flowfield(2).lambda_T = tangent_vel/(rotorsystem.rotor(2).omega*rotorsystem.rotor(2).R)*ones(res_psi,res_r); %%normalizing free stream tangential velocity by tip velocity - FOR LATER
    %% Sanity check - will have to do try/catch when choosing a value of trim

    if strcmpi(method,'leishman') && rotorsystem.rotor(2).pitch_root>collective_l && (rotorsystem.rotor(2).pitch_root>collective_l && strcmpi(rotorsystem.rotor(2).twist_type,'ideal')) || (rotorsystem.rotor(2).twistdeg>collective_l && strcmpi(rotorsystem.rotor(2).twist_type,'linear'))
        error('Rotor lower: Your twist is larger than your collective! Negative pitch angles near the tip will break Leishman')
    end
    
    %% Calculate inflow for bottom rotor
    
    separation = rotorsystem.params.interrotor_spacing*rotorsystem.rotor(1).R;
    
    %radial_contraction = getContraction(separation,CT) %Landgrebe
    
    radial_contraction = rotorsystem.params.rd;
    
    mean_downwash = mean(mean(lambda_i_u))*rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R; 
    %check if in hover to avoid singularities when lambda_i =0
    if axial_vel==0 && tangent_vel ==0
        skew_angle = 0; %for any lambda_i in hover
    else
        skew_angle = atan(tangent_vel/(axial_vel+mean_downwash));

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
    
    
    [Thrust_l, Torque_l, Power_l, CT_l, CP_l, dCT_l, dCP_l, lambda_l, Re_l, alpha_l,alpha_negatives_l, phi_l, ...
    F_l, weighted_swirl_ratio_l, FOM_l, velocity_dimensional_l,pitchdeg_l,r,dr,psi_l] = spinRotor(rotorsystem.rotor(2),atm,'CW',collective_l,lambda_P_bottom,flowfield(2).lambda_T,method,epsilon);

    %% Output Vars

        
    Thrust = Thrust_u+Thrust_l;
    Torque = Torque_u+Torque_l;
    Power = Power_u+Power_l;
    CT = CT_u+CT_l;
    CP = CP_u+CP_l;
    
    FOM_coax = rotorsystem.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper
    
    net_torque_coeff = (CP_u-CP_l)/CP;
    net_torque_dimensional = Torque_u-Torque_l;
    
    
    
end


    
 %% Display output text

if verbose
    disp('                                                  ')
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' rotor specs:'))
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' power coefficient [-] ',num2str(CP)))
    disp(strcat(rotorsystem.type,' thrust coefficient [-] ',num2str(CT)))
    disp(['Total thrust [N]',' ',num2str(Thrust)])
    disp(['Total power [W]',' ',num2str(Power)])
    disp(['Total torque [Nm]',' ',num2str(Torque)])
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
    diskPlot(r,psi_u,alpha_u,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
    subplot(2,3,2)
    diskPlot(r,psi_u,dCT_u,rotorsystem.state,'$dC_T$ [-]')
    subplot(2,3,3)
    diskPlot(r,psi_u,dCP_u,rotorsystem.state,'$dC_P$ [-]')
    subplot(2,3,4)
    diskPlot(r,psi_u,lambda_u,rotorsystem.state, '$\lambda_u$ [-]')
    subplot(2,3,5)
    diskPlot(r,psi_u,F_u,rotorsystem.state,'Prandtl [-]')
    subplot(2,3,6)
    diskPlot(r,psi_u,rad2deg(phi_u),rotorsystem.state,'$\phi$ [deg]')
    if strcmpi(rotorsystem.type,"coaxial")
        figure(100); clf;
        title('Lower Rotor')
        subplot(2, 3, 1)
        hold on
        diskPlot(r,psi_u,alpha_l,rotorsystem.state, '$\alpha$ [deg]') %add optional arguments
        subplot(2,3,2)
        diskPlot(r,psi_u,dCT_l,rotorsystem.state,'$dC_T$ [-]')
        subplot(2,3,3)
        diskPlot(r,psi_u,dCP_l,rotorsystem.state,'$dC_P$ [-]')
        subplot(2,3,4)
        diskPlot(r,psi_u,lambda_l,rotorsystem.state,'$\lambda_l$ [-]')
        subplot(2,3,5)
        diskPlot(r,psi_u,F_l,rotorsystem.state,'Prandtl [-]')
        subplot(2,3,6)
        diskPlot(r,psi_u,rad2deg(phi_l),rotorsystem.state,'$\phi$ [-]')
        figure(101)
        diskPlot(r,psi_u,lambda_P_bottom,rotorsystem.state,'Non-dimensional Inflow bottom rotor') 
    end    
    %% Validation plots
    if axial_vel==0 && tangent_vel==0
        if abs(CT-0.004)<0.00001 && net_torque_dimensional<0.1 && strcmpi(rotorsystem.name,"Harrington1")  
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






