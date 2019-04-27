
%maybe make this a function, like fly(collective1,collective2). This
%allows to trim the helicopter in hover. Whenever you change pitch and
%therefore thrust and torque it would be good to start the convergence of
%Fcf and lambda (in a separate function) with the values previously found.
init 

dr = 0.001;
r = dr:dr:1-5*dr; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

plot_inflow_loss = true;
plot_coeffs = true;

axial_vel = 0; %m/s
flowfield(1).lambda_inf = axial_vel/(rotor(1).rpm*2*pi*rotor(1).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
flowfield(2).lambda_inf = axial_vel/(rotor(2).rpm*2*pi*rotor(2).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity

epsilon = 0.001; %convergence accuracy

trim = 1; %1 means that both rotors have the same geometrical pitch, so same collective setting
pitchdeg = 14; %deg
rotor(1).pitch = deg2rad(pitchdeg)*ones(1,length(r)); %rad
rotor(2).pitch = trim*rotor(1).pitch; %rad

%% Converge Fcf and inflow ratio

[Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield(flowfield, r, epsilon, rotor); %this is now converging values for lambda (lambda_i+lambda_inf) and Fcf


%% Calculate thrust and torque coefficients

[CP, CT_u, CP_u, CT_l, CP_l, spanwise_coeffs] = get_coeffs(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, rotor,params, flowfield);

FOM_u = CT_u^(3/2)/(sqrt(2)*CP_u); %treated as a single rotor;
FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;

CT = CT_u + CT_l;

FOM_coax = params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper

T = atm.rho*pi*rotor(1).R^4*rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
P = atm.rho*pi*rotor(1).R^5*rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius

disp(['Net torque coefficient',' ',num2str((CP_u-CP_l)/CP_l)])
disp(['Coaxial system power coefficient',' ',num2str(CP)])
disp(['Coaxial system thrust coefficient',' ',num2str(CT)])
disp(['Total thrust [N]',' ',num2str(T)])
disp(['Total power [W]',' ',num2str(P)])

%% Plotting

if plot_inflow_loss
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
    
end

if plot_coeffs
    
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
    title('dCtl vs radius - Bottom')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCtl/dr')

    subplot(2, 2, 4)
    plot(r, spanwise_coeffs.dCp_l, 'b-.')
    title('dCpl vs radius - Bottom')
    xlabel('r/R')
    ylabel('Non-dimensionalized dCpl/dr')
end

