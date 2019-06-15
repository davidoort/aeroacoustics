function [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT_iter(rotorsystem,atm,epsilon,plots,verbose)
%{
AXIAL & FORWARD Flight
This function can calculate the spanwise thrust & power coefficients and
inflow ratio which allows to estimate the performance of a coaxial rotor in
axial or hovering flight. This function is run from BEMT.m when
tangential_vel = 0. It complements the BEMT_axial function in that it
avoids some of its assumptions by instead having a more iterative approach.

Inputs:
    rotorsystem - (struct object) with operational (state) and geometric
    variables of the coaxial or single rotor

    atm - (struct object) with atmospheric parameters such as air density

    epsilon - (scalar) convergence accuracy for CT-F-lambda iteration

    plots - (boolean) indicates if plots should be returned or not

    verbose - (boolean) indicates if output text should be displayed or not

Outputs:
    Thrust - (scalar) Total thrust of the coaxial rotor in AXIAL Flight

    Torque - (scalar) Total torque on the coaxial rotor in AXIAL Flight

    Power - (scalar) Total power consumed by the coaxial rotor in AXIAL
    Flight

    CT - (scalar) Thrust coefficient of coaxial (or single) rotor system

    CP - (scalar) Power/Torque coefficient of coaxial (or single) rotor
    system

    net_torque_coeff - (scalar) normalized difference in torque coefficient 
    between upper and lower rotor (0 for single rotor). Needed for trimming
    
    Plots (optionally)

Other m-files required: 

    covergeflowfield
    get_coeffs
    getChord
    getReynolds

MAT-files required: none

Literature referenced: 

    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (6) & (7) WITH SIGNIFICANT MODIFICATIONS.

    ETH slides! (momentum theory)

Assumptions:
    
    Same assumptions as BEMT_axial (the ones used in the paper by Leishman)
    with the exception of:

        - No small angle assumptions (for inflow angle)
        - Airfoil data lookup tables used (more accurate lift and drag 
        polars, possibly dependent on Reynolds number of blade section)

Ideas:

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
May 2019; Last revision: 10-June-2019
%}
%------------- BEGIN CODE --------------

debug = false;
%% Init
rmax = 0.99; %The dr and 0.99 is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.
if debug
    dr = 0.1;
    r_vec = rotorsystem.rotor(1).hub_radial_fraction+dr:dr:rmax; %non-dimensionalized by tip radius. Rotors have the same radius.
    dpsi_vec = linspace(0,2*pi,10); %length(r_vec)+1
else
    dr = 0.01;
    r_vec = rotorsystem.rotor(1).hub_radial_fraction+dr:dr:rmax; %non-dimensionalized by tip radius. Rotors have the same radius.
    dpsi_vec = linspace(0,2*pi,length(r_vec)+1); %length(r_vec)+1
end
dpsi = dpsi_vec(2)-dpsi_vec(1);


chord_vec = linspace(rotorsystem.rotor(1).root_chord,rotorsystem.rotor(1).tip_chord,length(r_vec));
%chord_vec(length(r_vec)+1:end) = []; %resizing of the array to match the length of r_vec. Couldn't find a more elegant way of doing this

r = repmat(r_vec,length(dpsi_vec),1);
chord = repmat(chord_vec,length(dpsi_vec),1);
psi = repmat(dpsi_vec',1,length(r_vec));
    
axial_vel = rotorsystem.state.axial_vel;
tangent_vel = rotorsystem.state.tangent_vel;

flowfield(1).lambda_P = axial_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(size(r)); %normalizing free stream axial velocity by tip velocity
flowfield(1).lambda_T = tangent_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(size(r)); %%normalizing free stream tangential velocity by tip velocity - FOR LATER

geometric_pitch = getPitch(r_vec,rotorsystem.rotor(1).twist_type,rotorsystem.rotor(1).twistdeg,rotorsystem.rotor(1).pitch_root);

pitchdeg = geometric_pitch+ones(size(geometric_pitch))*rotorsystem.state.collective;

pitchdeg = repmat(pitchdeg,length(dpsi_vec),1);

rotorsystem.rotor(1).pitch = deg2rad(pitchdeg); %rad - this might get more complicated when the function gets cyclic input. Or not

F_old = ones(size(r));
lambda_old = 0*ones(size(r));%0.01*ones(size(r));%flowfield(1).lambda_P; 
[phi,phi_old] = getInflowAngle(lambda_old,r,psi,flowfield(1).lambda_T);
dCTu_old = getdCT(rotorsystem.rotor(1),atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,flowfield(1).lambda_T);


%% Iterate upper rotor

err_old = 1;
while err_old>epsilon
    
    F = getPrandtlTipLoss(rotorsystem.rotor(1),phi_old,r);
    
    dCTu = getdCT(rotorsystem.rotor(1),atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,flowfield(1).lambda_T);
    
    lambda = getLambda(flowfield(1).lambda_P, dCTu, F_old, r, dr,dpsi);
                                   
    [phi_negatives,phi] = getInflowAngle(lambda,r,psi,flowfield(1).lambda_T);
    
    err = norm([F(~isnan(lambda))-F_old(~isnan(lambda)),lambda(~isnan(lambda))-lambda_old(~isnan(lambda)),dCTu(~isnan(lambda))-dCTu_old(~isnan(lambda)),phi(~isnan(lambda))-phi_old(~isnan(lambda))]);
    
    if abs(err-err_old)<1e-5
        warning('Error constant, stopping iteration')
        break
    end
    
    err_old = err;
    
    dCTu_old = dCTu;
    lambda_old = lambda;
    phi_old = phi;
    F_old = F;
end

%% Calculate nondimensional coefficients

dCPu = getdCP(rotorsystem.rotor(1),atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,flowfield(1).lambda_P,flowfield(1).lambda_T);

% Remove NaNs procedure
dCTu_sum = dCTu(~isnan(dCTu));
dCPu_sum = dCPu(~isnan(dCPu));

CT = sum(sum(dCTu_sum));
CP = sum(sum(dCPu_sum));




FOM = CT^(3/2)/(sqrt(2)*CP); %treated as a single rotor;

weighted_swirl_ratio = getSwirl(lambda,flowfield(1).lambda_P,flowfield(1).lambda_T,r,dr,dpsi,dCPu);

%% Reynolds - to be moved out of here

chord_eff = interp1(r(1,:),chord(1,:),0.7);
vel_eff = ((rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R*0.7).^2+(interp1(r(1,:),lambda(1,:),0.7)).^2).^(1/2);
Re_u_eff = getReynolds(vel_eff,chord_eff,atm.kin_visc);

v_tip = rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R;

velocity_dimensional = sqrt((v_tip*flowfield(1).lambda_T.*sin(psi)+v_tip*r).^2+lambda.^2);

Re_u =  getReynolds(velocity_dimensional,chord,atm.kin_visc); % to be used in the future in get2Dcoeffs

%% Calculate dimensional parameters

Thrust = atm.rho*pi*rotorsystem.rotor(1).R^4*rotorsystem.rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
Power = atm.rho*pi*rotorsystem.rotor(1).R^5*rotorsystem.rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
Torque = atm.rho*pi*rotorsystem.rotor(1).R^5*rotorsystem.rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius

net_torque_coeff = 0;

[alpha_u,alpha_neg] = getAoA(rotorsystem.rotor(1).pitch,phi);

%% Bottom rotor
if strcmpi(rotorsystem.type,"coaxial")
    %% Init
    error('Switch to single rotor!')
    CP_u = CP;
    CT_u = CT;
    FOM_u = FOM;
    
    flowfield(2).lambda_P = axial_vel/(coaxial.rotor(2).omega*rotorsystem.rotor(2).R)*ones(size(r)); %normalizing free stream axial velocity by tip velocity
    
    trim = coaxial.state.trim;
    coaxial.rotor(2).pitch = trim*coaxial.rotor(1).pitch; %rad
    

    %% Calculate nondimensional coeffs
    
    CT = CT_u+CT_l;
    
    FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;
    FOM_coax = coaxial.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper
    
    %% Calculate dimensional parameters
    Thrust = atm.rho*pi*coaxial.rotor(1).R^4*coaxial.rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
    Power = atm.rho*pi*coaxial.rotor(1).R^5*coaxial.rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
    Torque = atm.rho*pi*coaxial.rotor(1).R^5*coaxial.rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius
    
    net_torque_coeff = (CP_u-CP_l)/CP;
    
    net_torque_dimensional = net_torque_coeff*Torque;
    
    alpha_u = rad2deg(coaxial.rotor(1).pitch-atan((lambda_u+flowfield(1).lambda_inf)./r)); %this is being recomputed for the sake of not passing more arguments through get coeffs
    alpha_l = rad2deg(coaxial.rotor(2).pitch-atan(lambda_l+flowfield(2).lambda_inf)./r);
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
    end

end
        
%% Display output text

if verbose
    disp('                                                  ')
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' rotor specs:'))
    disp('--------------------------------------------------')
    disp(strcat(rotorsystem.type,' power coefficient [-] ',num2str(CP)))
    disp(strcat(rotorsystem.type,' thrust coefficient [-] ',num2str(CT)))
    disp(strcat(rotorsystem.type,' FOM [-] ',num2str(FOM)))
    disp(['Total thrust [N]',' ',num2str(Thrust)])
    disp(['Total power [W]',' ',num2str(Power)])
    disp(['Total torque [Nm]',' ',num2str(Torque)])
    disp(['Omega [rad/s]',' ',num2str(rotorsystem.rotor(1).omega)])
    disp(['Swirl upper rotor [rad/s] ',num2str(weighted_swirl_ratio*rotorsystem.rotor(1).omega)])
    disp(['Max/Min AoA upper rotor [deg] ',num2str(max(max(alpha_u))),' / ',num2str(min(min(alpha_u)))])
    disp(['Max/Min/0.7R Re number upper [-] ',num2str(max(max(Re_u))),' / ',num2str(min(min(Re_u))),' / ',num2str(Re_u_eff)])
    if strcmpi(rotorsystem.type,"coaxial")
        disp(['Net torque coefficient (u-l)/l [-]',' ',num2str(net_torque_coeff)])
        disp(['Net torque (u-l) [Nm]',' ',num2str(net_torque_dimensional)])
        disp(['Pitch upper/lower rotor [deg] ',num2str(pitchdeg),' / ',num2str(rad2deg(coaxial.rotor(2).pitch(1)))])
        disp(['Max/Min AoA lower rotor [deg] ',num2str(max(alpha_l)),' / ',num2str(min(alpha_l))])
        disp(['Max/Min/0.7R Re number lower [-] ',num2str(max(Re_l)),' / ',num2str(min(Re_l)),' / ',num2str(Re_l_eff)])
    else
        disp(['Pitch rotor [deg] ',num2str(pitchdeg(1,:))])        
    end
    
end

%% Plots

if plots
    figure(1); clf;
    subplot(2, 3, 1)
    hold on
    diskPlot(r,psi,alpha_u,'alpha', rotorsystem.rotor(1).hub_radial_fraction) %add optional arguments
    subplot(2,3,2)
    diskPlot(r,psi,dCTu,'dCT', rotorsystem.rotor(1).hub_radial_fraction)
    subplot(2,3,3)
    diskPlot(r,psi,dCPu,'dCP', rotorsystem.rotor(1).hub_radial_fraction)
    subplot(2,3,4)
    diskPlot(r,psi,lambda-flowfield(1).lambda_P,'induced inflow', rotorsystem.rotor(1).hub_radial_fraction)
    subplot(2,3,5)
    diskPlot(r,psi,F,'Prandtl', rotorsystem.rotor(1).hub_radial_fraction)
    subplot(2,3,6)
    diskPlot(r,psi,rad2deg(phi),'Phi', rotorsystem.rotor(1).hub_radial_fraction)
end

end

