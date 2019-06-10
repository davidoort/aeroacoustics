function [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT_axial(rotor,atm,epsilon,plots,verbose)
%{
AXIAL Flight
This function can calculate the spanwise thrust & power coefficients and
inflow ratio which allows to estimate the performance of a coaxial rotor in
axial or hovering flight. This function is run from BEMT.m when
tangential_vel = 0.

Inputs:
    rotor - (struct object) with operational (state) and geometric
    variables of the coaxial or single rotor

    atm - (struct object) with atmospheric parameters such as air density

    epsilon - (scalar) convergence accuracy for F-lambda iteration

    plots - (boolean) indicates if plots should be returned or not

Outputs:
    Thrust - (scalar) Total thrust of the coaxial rotor in AXIAL Flight

    Torque - (scalar) Total torque on the coaxial rotor in AXIAL Flight

    Power - (scalar) Total power consumed by the coaxial rotor in AXIAL
    Flight

    Plots (optionally)

Other m-files required: 

    covergeflowfield
    get_coeffs

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (6) & (7)

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions:
    In the individual functions


Ideas:
    %maybe make this a function, like fly(collective1,collective2). This
    %allows to trim the helicopter in hover. Whenever you change pitch and
    %therefore thrust and torque it would be good to start the convergence of
    %Fcf and lambda (in a separate function) with the values previously found.


Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
May 2019; Last revision: 9-June-2019
%}
%------------- BEGIN CODE --------------

dr = 0.01;
r = rotor.rotor(1).hub_radial_fraction+dr:dr:0.99; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

if strcmpi(rotor.type,"single")
    single = rotor;
    axial_vel = single.state.axial_vel;
    
    
    flowfield(1).lambda_inf = axial_vel/(single.rotor(1).omega*single.rotor(1).R)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
    
    pitchdeg = single.state.pitchdeg;
    
    single.rotor(1).pitch = deg2rad(pitchdeg)*ones(1,length(r)); %rad
    
    %% Converge Fcf and inflow ratio
    
    [Fcf_u, lambda_u] = convergeflowfield_single(flowfield, r, epsilon, single.rotor); %this is now converging values for lambda (lambda_i+lambda_inf) and Fcf
    
    %% Calculate thrust and torque coefficients
    
    [CT, CP, spanwise_coeffs, weighted_swirl_ratio] = get_coeffs_single(Fcf_u, lambda_u, r, dr, single, flowfield);
    
    
    Thrust = atm.rho*pi*single.rotor(1).R^4*single.rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
    Power = atm.rho*pi*single.rotor(1).R^5*single.rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
    Torque = atm.rho*pi*single.rotor(1).R^5*single.rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius
    
    net_torque_coeff = 0;
    
    alpha_u = rad2deg(single.rotor(1).pitch-atan((lambda_u+flowfield(1).lambda_inf)./r)); 
    
    if verbose
        disp(['Single rotor power coefficient [-]',' ',num2str(CP)])
        disp(['Single rotor thrust coefficient [-]',' ',num2str(CT)])
        disp(['Total thrust [N]',' ',num2str(Thrust)])
        disp(['Total power [W]',' ',num2str(Power)])
        disp(['Total torque [Nm]',' ',num2str(Torque)])
        disp(['Omega [rad/s]',' ',num2str(single.rotor(1).omega)])
        disp(['Swirl rotor [rad/s] ',num2str(weighted_swirl_ratio*single.rotor(1).omega)])
        disp(['Pitch rotor [deg] ',num2str(pitchdeg)])
        disp(['Max/Min AoA upper rotor [deg] ',num2str(max(alpha_u)),' / ',num2str(min(alpha_u))])
    end
    
else
    
    coaxial = rotor;
        
    axial_vel = coaxial.state.axial_vel;
    
    flowfield(1).lambda_inf = axial_vel/(coaxial.rotor(1).omega*coaxial.rotor(1).R)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
    flowfield(2).lambda_inf = axial_vel/(coaxial.rotor(2).omega*coaxial.rotor(2).R)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
    
    pitchdeg = coaxial.state.pitchdeg;
    
    trim = coaxial.state.trim;
    
    
    coaxial.rotor(1).pitch = deg2rad(pitchdeg)*ones(1,length(r)); %rad
    coaxial.rotor(2).pitch = trim*coaxial.rotor(1).pitch; %rad
    
    %% Converge Fcf and inflow ratio
    
    [Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield(flowfield, r, epsilon, coaxial); %this is now converging values for lambda (lambda_i+lambda_inf) and Fcf
    
    
    %% Calculate thrust and torque coefficients
    
    [CP, CT_u, CP_u, CT_l, CP_l, spanwise_coeffs] = get_coeffs(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, coaxial, flowfield);
    
    CT = CT_u + CT_l;
    
    %% Calculate FOM
    FOM_u = CT_u^(3/2)/(sqrt(2)*CP_u); %treated as a single rotor;
    FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;
    FOM_coax = coaxial.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper
    
    %% Calculate dimensional things
    Thrust = atm.rho*pi*coaxial.rotor(1).R^4*coaxial.rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
    Power = atm.rho*pi*coaxial.rotor(1).R^5*coaxial.rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
    Torque = atm.rho*pi*coaxial.rotor(1).R^5*coaxial.rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius
    
    net_torque_coeff = (CP_u-CP_l)/CP;
    
    net_torque_dimensional = net_torque_coeff*Torque;
    
    alpha_u = rad2deg(coaxial.rotor(1).pitch-atan((lambda_u+flowfield(1).lambda_inf)./r)); %this is being recomputed for the sake of not passing more arguments through get coeffs
    alpha_l = rad2deg(coaxial.rotor(2).pitch-atan(lambda_l+flowfield(2).lambda_inf)./r);
    
    %% Reynolds number calculation
    
    omega_u = coaxial.rotor(1).omega;
    omega_l = coaxial.rotor(2).omega;
    
    R_u = coaxial.rotor(1).R;
    R_l = coaxial.rotor(2).R;
    
    
    chord_u         = getChord(r,coaxial.rotor(1).hub_radial_fraction,coaxial.rotor(1).tip_chord,coaxial.rotor(1).root_chord);
    chord_u_eff     = getChord(0.7,coaxial.rotor(1).hub_radial_fraction,coaxial.rotor(1).tip_chord,coaxial.rotor(1).root_chord);
    chord_l         = getChord(r,coaxial.rotor(2).hub_radial_fraction,coaxial.rotor(2).tip_chord,coaxial.rotor(2).root_chord);
    chord_l_eff     = getChord(0.7,coaxial.rotor(2).hub_radial_fraction,coaxial.rotor(2).tip_chord,coaxial.rotor(2).root_chord);
    
    vel_u     = ((omega_u*R_u*r).^2+(lambda_u+flowfield(1).lambda_inf).^2).^(1/2);
    vel_u_eff = ((omega_u*R_u*0.7).^2+(lambda_u(abs(r - 0.7)<=dr)+flowfield(1).lambda_inf(abs(r - 0.7)<=dr)).^2).^(1/2);
    vel_l     = ((omega_l*R_l*r).^2+(lambda_l+flowfield(2).lambda_inf).^2).^(1/2);
    vel_l_eff = ((omega_l*R_l*0.7).^2+(lambda_l(abs(r - 0.7)<=dr)+flowfield(2).lambda_inf(abs(r - 0.7)<=dr)).^2).^(1/2);
    
    Re_u = getReynolds(vel_u,chord_u,atm.kin_visc);
    Re_u_eff = mean(getReynolds(vel_u_eff,chord_u_eff,atm.kin_visc));
    Re_l = getReynolds(vel_l,chord_l,atm.kin_visc);
    Re_l_eff = mean(getReynolds(vel_l_eff,chord_l_eff,atm.kin_visc));
    
    %% Display output text
    if verbose
        disp(['Net torque coefficient (u-l)/l [-]',' ',num2str(net_torque_coeff)])
        disp(['Net torque (u-l) [Nm]',' ',num2str(net_torque_dimensional)])
        disp(['Coaxial system power coefficient [-]',' ',num2str(CP)])
        disp(['Coaxial system thrust coefficient [-]',' ',num2str(CT)])
        disp(['Total thrust [N]',' ',num2str(Thrust)])
        disp(['Total power [W]',' ',num2str(Power)])
        disp(['Pitch upper/lower rotor [deg] ',num2str(pitchdeg),' / ',num2str(rad2deg(coaxial.rotor(2).pitch(1)))])
        disp(['Max/Min AoA upper rotor [deg] ',num2str(max(alpha_u)),' / ',num2str(min(alpha_u))])
        disp(['Max/Min AoA lower rotor [deg] ',num2str(max(alpha_l)),' / ',num2str(min(alpha_l))])
        disp(['Max/Min/0.7R Re number upper [-] ',num2str(max(Re_u)),' / ',num2str(min(Re_u)),' / ',num2str(Re_u_eff)])
        disp(['Max/Min/0.7R Re number lower [-] ',num2str(max(Re_l)),' / ',num2str(min(Re_l)),' / ',num2str(Re_l_eff)])
    end
    %% Plotting
    
    if abs(CT-0.004)<0.00001 && net_torque_dimensional<0.1 && strcmpi(coaxial.name,"Harrington1")  && strcmpi(coaxial.type,"coaxial")
       
       data_inflow_upper = readmatrix('H1_inflow_FVM_fig10a.csv');
       r1 = data_inflow_upper(:,1); lambda1 = data_inflow_upper(:,2);
       data_inflow_lower = readmatrix('H1_inflow_FVM_fig10b.csv'); 
       r2 = data_inflow_lower(:,1); lambda2 = data_inflow_lower(:,2);
       
       data_dCT_upper = readmatrix('H1_dCT_FVM_fig11a.csv'); 
       r3 = data_dCT_upper(:,1); dCT1 = data_dCT_upper(:,2);
       data_dCT_lower = readmatrix('H1_dCT_FVM_fig11b.csv'); 
       r4 = data_dCT_lower(:,1); dCT2 = data_dCT_lower(:,2);
       
       data_dCP_upper = readmatrix('H1_dCP_FVM_fig12a.csv');
       r5 = data_dCP_upper(:,1); dCP1 = data_dCP_upper(:,2);    
       data_dCP_lower = readmatrix('H1_dCP_FVM_fig12b.csv'); 
       r6 = data_dCP_lower(:,1); dCP2 = data_dCP_lower(:,2); 
       
    else
        
        r1=[];r2=[];r3=[];r4=[];r5=[];r6=[];
        lambda1=[];lambda2=[];dCT1=[];dCT2=[];dCP1=[];dCP2=[];
        
    end
    
    
    if plots
        
        %%%%%
        
        figure(1); clf;
        subplot(2, 3, 1)
        hold on
        scatter(r1,lambda1)
        plot(r, lambda_u, 'b-.')
        title('Inflow ratio vs radius - Top')
        xlabel('r/R')
        ylabel('Inflow ratio')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        subplot(2, 3, 2)
        plot(r, Fcf_u, 'b-.')
        title('Prandtl tip loss vs radius - Top')
        xlabel('r/R')
        ylabel('Prandtl tip loss')
         
        subplot(2, 3, 3)
        plot(r, alpha_u, 'b-.')
        title('AoA vs radius - Top')
        xlabel('r/R')
        ylabel('$\alpha$ [deg]','interpreter','latex')
        ylim([-20 Inf])
        
        subplot(2, 3, 4)
        hold on
        scatter(r2,lambda2)
        plot(r, lambda_l, 'b-.')
        title('Inflow ratio vs radius - Bottom')
        xlabel('r/R')
        ylabel('Inflow ratio')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        subplot(2, 3, 5)
        plot(r, Fcf_l, 'b-.')
        title('Prandtl tip loss vs radius - Bottom')
        xlabel('r/R')
        ylabel('Prandtl tip loss')

        subplot(2, 3, 6)
        plot(r, alpha_l, 'b-.')
        title('AoA vs radius - Bottom')
        xlabel('r/R')
        ylabel('$\alpha$ [deg]','interpreter','latex')
        ylim([-20 Inf])
        
        %%%%%%%
        
        figure(2); clf;
        subplot(2, 2, 1)
        hold on
        scatter(r3,dCT1)
        plot(r, spanwise_coeffs.dCt_u, 'b-.')
        title('dCTu vs radius - Top')
        xlabel('r/R')
        ylabel('Non-dimensionalized dCTu/dr')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        subplot(2, 2, 2)
        hold on
        scatter(r5,dCP1)
        plot(r, spanwise_coeffs.dCp_u, 'b-.')
        title('dCpu vs radius - Top')
        xlabel('r/R')
        ylabel('Non-dimensionalized dCpu/dr')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        subplot(2, 2, 3)
        hold on
        scatter(r4,dCT2)
        plot(r, spanwise_coeffs.dCt_l, 'b-.')
        title('dCTl vs radius - Bottom')
        xlabel('r/R')
        ylabel('Non-dimensionalized dCtl/dr')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        subplot(2, 2, 4)
        hold on
        scatter(r6,dCP2)
        plot(r, spanwise_coeffs.dCp_l, 'b-.')
        title('dCpl vs radius - Bottom')
        xlabel('r/R')
        ylabel('Non-dimensionalized dCpl/dr')
        xlim([0.2 1])
        legend('FVM','BEMT')
        
        
        %% Poster plots
        
        
        % Disk-plot of velocity tangential velocity at the blade in forward flight
        
%         phi = linspace(0,2*pi);
%         radius = linspace(0.15,1);
%         
%         [PHI,R] = meshgrid(phi,radius);
%         
%         vel = sin(PHI)+R;
%         
%         figure(3)
%         surf(R*cos(PHI),R*sin(PHI),vel)
%         title('Velocity along $\vec{y_b}$','interpreter','latex')
%         %xlabel('r/R')
%         %ylabel('Non-dimensionalized dCpl/dr')
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         set(gca,'visible','off')
%         colorbar
%         
%         % Combine previous plots to make them denser (ideally they would have been combined with the plots in Leishman's paper)
%         
%         %set(gca,'fontsize',30)
%         
%         figure(4); clf;
%         
%         hold on;
%         plot(r, lambda_u, '-.', 'LineWidth',2)
%         plot(r, lambda_l, '-.', 'LineWidth',2)
%         lgd=legend('Upper rotor','Lower rotor');
%         lgd.FontSize = 14;
%         %title('Inflow ratio vs radius')
%         xlabel('r/R','FontSize',14)
%         ylabel('Inflow ratio','FontSize',14)
%         
%         figure(5); clf;
%         hold on;
%         plot(r, spanwise_coeffs.dCt_u, '-.', 'LineWidth',2)
%         plot(r, spanwise_coeffs.dCt_l, '-.', 'LineWidth',2)
%         lgd=legend('Upper rotor','Lower rotor');
%         lgd.FontSize = 14;
%         %title('dCTu vs radius')
%         xlabel('r/R','FontSize',14)
%         ylabel('Non-dimensionalized dCT/dr','FontSize',14)
%         
%         figure(6); clf;
%         hold on;
%         plot(r, spanwise_coeffs.dCp_u, '-.', 'LineWidth',2)
%         plot(r, spanwise_coeffs.dCp_l, '-.', 'LineWidth',2)
%         lgd=legend('Upper rotor','Lower rotor');
%         lgd.FontSize = 14;
%         %title('dCTu vs radius')
%         xlabel('r/R','FontSize',14)
%         ylabel('Non-dimensionalized dCP/dr','FontSize',14)
        
        
    end
    
end

end

