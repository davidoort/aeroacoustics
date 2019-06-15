function [Thrust, Torque, Power, CT, CP, net_torque_coeff] = BEMT(rotorsystem,atm,epsilon,plots,verbose,method,debug)


%init is done inside the spinRotor function, but you have to specify which
%rotor (geometry) has to be spun and in which direction it has to be spun, with what
%collective (and cyclic later), lambda_P, lambda_T, method (leishman or airfoil)

if debug
    res_r = 6;
    res_psi = 10;
else
    res_r = 100;
    res_psi = 130; %careful of putting them equal
end

axial_vel = rotorsystem.state.axial_vel;
tangent_vel = rotorsystem.state.tangent_vel;

flowfield(1).lambda_P = axial_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %normalizing free stream axial velocity by tip velocity
flowfield(1).lambda_T = tangent_vel/(rotorsystem.rotor(1).omega*rotorsystem.rotor(1).R)*ones(res_psi,res_r); %%normalizing free stream tangential velocity by tip velocity - FOR LATER

collective_u = rotorsystem.state.collective;

[Thrust, Torque, Power, CT, CP, dCT, dCP, lambda, Re_u, alpha_u, phi, ...
    F, weighted_swirl_ratio, FOM, velocity_dimensional,pitchdeg,r,psi] = spinRotor(rotorsystem.rotor(1),atm,'CCW',collective_u,flowfield(1).lambda_P,flowfield(1).lambda_T,method,epsilon,debug);

net_torque_coeff = 0;


if strcmpi(rotorsystem.type,"coaxial")
    error('Switch to single rotor!')
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
    disp(['Max/Min/0.7R Re number upper [-] ',num2str(max(max(Re_u))),' / ',num2str(min(min(Re_u))),' / ',num2str(interp1(r(1,:),Re_u(1,:),0.7))])
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

%DISK PLOTS

if plots
    figure(1); clf;
    subplot(2, 3, 1)
    hold on
    diskPlot(r,psi,alpha_u,'alpha') %add optional arguments
    subplot(2,3,2)
    diskPlot(r,psi,dCT,'dCT')
    subplot(2,3,3)
    diskPlot(r,psi,dCP,'dCP')
    subplot(2,3,4)
    diskPlot(r,psi,lambda-flowfield(1).lambda_P,'induced inflow')
    subplot(2,3,5)
    diskPlot(r,psi,F,'Prandtl')
    subplot(2,3,6)
    diskPlot(r,psi,rad2deg(phi),'Phi')
end



end

%For later
%{






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
    plot(r, AoA, 'b-.')
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

%}





