function [Thrust, Torque, Power] = BEMT_FF(coaxial,atm,epsilon,plots)
%{
FORWARD Flight
This function can calculate the spanwise thrust & power coefficients and
inflow ratio which allows to estimate the performance of a coaxial rotor in
forward, axial or hovering flight (by discretizing into disk elements). 
This function is run from BEMT.m when there is a tangential velocity
component.

Inputs:
    coaxial - (struct object) with operational (state) and geometric
    variables of the coaxial rotor

    atm - (struct object) with atmospheric parameters such as air density

    epsilon - (scalar) convergence accuracy for F-lambda iteration

    plots - (boolean) indicates if plots should be returned or not

Outputs:
    Thrust - (scalar) Total thrust of the coaxial rotor in FF

    Torque - (scalar) Total torque on the coaxial rotor in FF

    Power - (scalar) Total power consumed by the coaxial rotor in FF
    

    Plots (optionally)

Other m-files required: 

    covergeflowfield_FF
    get_coeffs_FF

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

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
May 2019; Last revision: 26-May-2019
%}



%------------- BEGIN CODE --------------


%% Discretization in r and phi

dr = 0.001;
r = dr:dr:1-5*dr; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

dphi = 0.01;
psi_arr = 0:dphi:2*pi;

%% Converge upper rotor (this could be inside a big function like converge flowfield)

for psi = psi_arr
    %% Define spanwise variables (could be a function of phi)
    rotor = coaxial.rotor;

    axial_vel = coaxial.state.axial_vel;

    flowfield(1).lambda_P = axial_vel/(rotor(1).rpm*2*pi*rotor(1).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
   

    pitchdeg = coaxial.state.pitchdeg;

    trim = coaxial.state.trim;

    coaxial.rotor(1).pitch = deg2rad(pitchdeg)*ones(1,length(r)); %rad - if you had cyclic input and twist (Talaria implementation) this will have to be changed
    
    
    pitch_u = rotor(1).pitch;

    
    %% Iteration initialization
    
    Fcf0_u = zeros(1,length(r)); %dummy, to start iteration
    Fcf_u = ones(1,length(r)); %"Fcf started with an initial value of 1" - in standard procedure also by Leishman
   
    lambda0_u = zeros(1,length(r)); %dummy, to start iteration
    
    lambda_tot_u = get_lambda_up_FF(Fcf_u,r,psi, pitch_u,rotor,flowfield); %remember that this is lambda_u_induced+lambda_P
    
    
    %% converge upper rotor - can be done without converging bottom rotor (assumption of no influence of bottom on top)
    i = 0;
    while norm(Fcf_u-Fcf0_u)>epsilon || norm(lambda_tot_u-lambda0_u)>epsilon
        Fcf0_u = Fcf_u;
        lambda0_u = lambda_tot_u;
        
        Fcf_u = Prandtl_tip_loss_FF(r,lambda_tot_u,rotor(1));
        lambda_tot_u = get_lambda_up_FF(Fcf_u,r,pitch_u,rotor,flowfield);
        i = i+1;
    end
    
    
    
end


%% Converge lower rotor

flowfield(2).lambda_P = axial_vel/(rotor(2).rpm*2*pi*rotor(2).R/60)*ones(1,length(r)); %normalizing free stream axial velocity by tip velocity
coaxial.rotor(2).pitch = trim*coaxial.rotor(1).pitch; %rad - for the Talaria implementation this will work in a similar way since the cyclic is coupled and you only have freedom in the collective



pitch_l = rotor(2).pitch;
Fcf0_l = zeros(length(r),length(psi)); %dummy, to start iteration
Fcf_l = ones(length(r),length(psi)); %"Fcf started with an initial value of 1"





%% Old

% Converge Fcf and inflow ratio

%[Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield_FF(flowfield, r, phi, epsilon, coaxial); %this is now converging values for lambda (lambda_i+lambda_P) and Fcf




% %% Calculate thrust and torque coefficients
% 
% [CP, CT_u, CP_u, CT_l, CP_l, spanwise_coeffs] = get_coeffs_FF(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, coaxial, flowfield);
% 
% FOM_u = CT_u^(3/2)/(sqrt(2)*CP_u); %treated as a single rotor;
% FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;
% 
% CT = CT_u + CT_l;
% 
% FOM_coax = coaxial.params.kappaint*(CT_u^(3/2)+CT_l^(3/2))/(sqrt(2)*(CP_u+CP_l)); %from robust control paper
% 
% Thrust = atm.rho*pi*rotor(1).R^4*rotor(1).omega^2*CT; %Using rotor 1 radius=rotor2 radius
% Power = atm.rho*pi*rotor(1).R^5*rotor(1).omega^3*CP; %Using rotor 1 radius=rotor2 radius
% Torque = atm.rho*pi*rotor(1).R^5*rotor(1).omega^2*CP; %Using rotor 1 radius=rotor2 radius
% 
% disp(['Net torque coefficient',' ',num2str((CP_u-CP_l)/CP_l)])
% disp(['Coaxial system power coefficient',' ',num2str(CP)])
% disp(['Coaxial system thrust coefficient',' ',num2str(CT)])
% disp(['Total thrust [N]',' ',num2str(Thrust)])
% disp(['Total power [W]',' ',num2str(Power)])


end

