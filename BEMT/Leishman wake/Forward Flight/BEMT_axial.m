function [Thrust, Torque, Power] = BEMT_axial(coaxial,atm,epsilon,plots)
%maybe make this a function, like fly(collective1,collective2). This
%allows to trim the helicopter in hover. Whenever you change pitch and
%therefore thrust and torque it would be good to start the convergence of
%Fcf and lambda (in a separate function) with the values previously found.

dr = 0.001;
r = dr:dr:1-5*dr; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

rotor = coaxial.rotor;

axial_vel = coaxial.state.axial_vel;

flowfield(1).lambda_inf = axial_vel/(rotor(1).rpm*2*pi*rotor(1).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
flowfield(2).lambda_inf = axial_vel/(rotor(2).rpm*2*pi*rotor(2).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity

pitchdeg = coaxial.state.pitchdeg;

trim = coaxial.state.trim;


coaxial.rotor(1).pitch = deg2rad(pitchdeg)*ones(1,length(r)); %rad
coaxial.rotor(2).pitch = trim*coaxial.rotor(1).pitch; %rad

%% Converge Fcf and inflow ratio

[Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield(flowfield, r, epsilon, coaxial); %this is now converging values for lambda (lambda_i+lambda_inf) and Fcf


%% Calculate thrust and torque coefficients

[CP, CT_u, CP_u, CT_l, CP_l, spanwise_coeffs] = get_coeffs(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, coaxial, flowfield);

FOM_u = CT_u^(3/2)/(sqrt(2)*CP_u); %treated as a single rotor;
FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;

CT = CT_u + CT_l;

FOM_coax = coaxial.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper

Thrust = atm.rho*pi*rotor(1).R^4*rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
Power = atm.rho*pi*rotor(1).R^5*rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
Torque = atm.rho*pi*rotor(1).R^5*rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius

disp(['Net torque coefficient',' ',num2str((CP_u-CP_l)/CP_l)])
disp(['Coaxial system power coefficient',' ',num2str(CP)])
disp(['Coaxial system thrust coefficient',' ',num2str(CT)])
disp(['Total thrust [N]',' ',num2str(Thrust)])
disp(['Total power [W]',' ',num2str(Power)])

%% Plotting

if plots
    figure(1); clf;
    subplot(2, 2, 1)
    plot(r, lambda_u, 'b-.')
    title('Inflow ratio vs radius - Top')
    xlabel('r/R')
    ylabel('Inflow ratio')

    subplot(2, 2, 2)
    plot(r, Fcf_u, 'b-.')
    title('Prandtl tip loss vs radius - Top')
    xlabel('r/R')
    ylabel('Prandtl tip loss')

    subplot(2, 2, 3)
    plot(r, lambda_l, 'b-.')
    title('Inflow ratio vs radius - Bottom')
    xlabel('r/R')
    ylabel('Inflow ratio')

    subplot(2, 2, 4)
    plot(r, Fcf_l, 'b-.')
    title('Prandtl tip loss vs radius - Bottom')
    xlabel('r/R')
    ylabel('Prandtl tip loss')
    
    figure(2); clf;
    subplot(2, 2, 1)
    plot(r, spanwise_coeffs.dCt_u, 'b-.')
    title('dCTu vs radius - Top')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCTu/dr')

    subplot(2, 2, 2)
    plot(r, spanwise_coeffs.dCp_u, 'b-.')
    title('dCpu vs radius - Top')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCpu/dr')

    subplot(2, 2, 3)
    plot(r, spanwise_coeffs.dCt_l, 'b-.')
    title('dCTl vs radius - Bottom')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCtl/dr')

    subplot(2, 2, 4)
    plot(r, spanwise_coeffs.dCp_l, 'b-.')
    title('dCpl vs radius - Bottom')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCpl/dr')
end


%% Poster plots 


% Disk-plot of velocity tangential velocity at the blade in forward flight

phi = linspace(0,2*pi);
radius = linspace(0.15,1);

[PHI,R] = meshgrid(phi,radius);

vel = sin(PHI)+R;

figure(3)
surf(R*cos(PHI),R*sin(PHI),vel)
title('Velocity along $\vec{y_b}$','interpreter','latex')
%xlabel('r/R')
%ylabel('Non-dimensionalized dCpl/dr')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'visible','off')
colorbar

% Combine previous plots to make them denser (ideally they would have been combined with the plots in Leishman's paper)

%set(gca,'fontsize',30)

figure(4); clf;

hold on;
plot(r, lambda_u, '-.', 'LineWidth',2)
plot(r, lambda_l, '-.', 'LineWidth',2)
lgd=legend('Upper rotor','Lower rotor');
lgd.FontSize = 14;
%title('Inflow ratio vs radius')
xlabel('r/R','FontSize',14)
ylabel('Inflow ratio','FontSize',14)

figure(5); clf;
hold on;
plot(r, spanwise_coeffs.dCt_u, '-.', 'LineWidth',2)
plot(r, spanwise_coeffs.dCt_l, '-.', 'LineWidth',2)
lgd=legend('Upper rotor','Lower rotor');
lgd.FontSize = 14;
%title('dCTu vs radius')
xlabel('r/R','FontSize',14)
ylabel('Non-dimensionalized dCT/dr','FontSize',14)

figure(6); clf;
hold on;
plot(r, spanwise_coeffs.dCp_u, '-.', 'LineWidth',2)
plot(r, spanwise_coeffs.dCp_l, '-.', 'LineWidth',2)
lgd=legend('Upper rotor','Lower rotor');
lgd.FontSize = 14;
%title('dCTu vs radius')
xlabel('r/R','FontSize',14)
ylabel('Non-dimensionalized dCP/dr','FontSize',14)






end

