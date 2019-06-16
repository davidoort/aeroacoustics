function [Thrust, Torque, Power, CT, CP, dCT, dCP, lambda, Re, AoA,alpha_negatives, phi, F, weighted_swirl_ratio,FOM_single,velocity_dimensional,pitchdeg,r,dr,psi] = spinRotor(rotor,atm,spin_dir,collective,lambda_P,lambda_T,method,epsilon)
%{
spinRotor
Modular function that can spin any rotor wit

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
June 2019; Last revision: 15-June-2019
%}
%------------- BEGIN CODE --------------


if strcmpi(spin_dir,'CCW')
    spin = 1;
elseif strcmpi(spin_dir,'CW')
    spin = -1;
end

%% Init
rmax = 0.99; %The dr and 0.99 is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

r_vec = linspace(rotor.hub_radial_fraction,rmax,size(lambda_P,2));
psi_vec = linspace(0,spin*2*pi,size(lambda_P,1));

dr = r_vec(2)-r_vec(1);
dpsi = spin*(psi_vec(2)-psi_vec(1));

chord_vec = linspace(rotor.root_chord,rotor.tip_chord,length(r_vec));
%chord_vec(length(r_vec)+1:end) = []; %resizing of the array to match the length of r_vec. Couldn't find a more elegant way of doing this

r = repmat(r_vec,length(psi_vec),1);
chord = repmat(chord_vec,length(psi_vec),1);
psi = repmat(psi_vec',1,length(r_vec));    

geometric_pitch = getPitch(r_vec,rotor.twist_type,rotor.twistdeg,rotor.pitch_root);

pitchdeg = geometric_pitch+ones(size(geometric_pitch))*collective;

pitchdeg = repmat(pitchdeg,length(psi_vec),1);

rotor.pitch = deg2rad(pitchdeg); %rad - this might get more complicated when the function gets cyclic input. Or not

%% Iteration
if strcmpi(method,'airfoil')
    F_old = ones(size(r));
    lambda_old = 0.0001*ones(size(r));%0.01*ones(size(r));%lambda_P;
    [phi,phi_old] = getInflowAngle(lambda_old,lambda_P,r,psi,lambda_T);
    dCTu_old = getdCT(rotor,atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,lambda_T);
    
    err_old = 1;
    while err_old>epsilon
        
        F = getPrandtlTipLoss(rotor,phi_old,r);
        
        dCT = getdCT(rotor,atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,lambda_T);
        
        lambda = getLambda(lambda_P, dCT, F_old, r, dr,dpsi);
        
        [phi_negative,phi] = getInflowAngle(lambda,lambda_P,r,psi,lambda_T);
        
        err = norm([F(~isnan(lambda))-F_old(~isnan(lambda)),lambda(~isnan(lambda))-lambda_old(~isnan(lambda)),dCT(~isnan(lambda))-dCTu_old(~isnan(lambda)),phi(~isnan(lambda))-phi_old(~isnan(lambda))]);
        
        if abs(err-err_old)<=1e-4
            warning('Error constant, stopping iteration')
            break
        end
        
        err_old = err;
        
        dCTu_old = dCT;
        lambda_old = lambda;
        phi_old = phi;
        F_old = F;
    end
    
    
    dCP = getdCP(rotor,atm,phi_old,r,dr,psi,dpsi,chord,lambda_old,lambda_P,lambda_T);

elseif strcmpi(method,'leishman')
    
    
    
    F_old = ones(size(r));
    lambda_old = getLambda_Leish(rotor,lambda_P,lambda_T,F_old,r,psi);
    [phi_negative,phi_old] = getInflowAngle(lambda_old,lambda_P,r,psi,lambda_T);
    
    err = 1;
    while err>epsilon
        
        F = getPrandtlTipLoss(rotor,phi_old,r);
        
        %lambda = getLambda_Leish(lambda_P, dCTu,F_old, r, dr,dpsi);
        
        lambda = getLambda_Leish(rotor,lambda_P,lambda_T,F,r,psi);
        
        [phi_negative,phi] = getInflowAngle(lambda,lambda_P,r,psi,lambda_T);
        
        err = norm([F-F_old,lambda-lambda_old,phi-phi_old]);
        
        lambda_old = lambda;
        phi_old = phi;
        F_old = F;
    end
    
    dCT = getdCT_Leish(rotor,F,lambda,lambda_P,lambda_T,phi,r,dr,psi,dpsi);
    
    dCP = getdCP_Leish(rotor,dCT,phi,r,dr,psi,dpsi,lambda,lambda_T);

    
end

%% Calculate nondimensional coefficients

% Remove NaNs procedure
dCTu_sum = dCT(~isnan(dCT));
dCPu_sum = dCP(~isnan(dCP));

CT = sum(sum(dCTu_sum));
CP = sum(sum(dCPu_sum));

weighted_swirl_ratio = getSwirl(lambda,lambda_P,lambda_T,r,dr,dpsi,dCP);

FOM_single = CT^(3/2)/(sqrt(2)*CP); %treated as a single rotor;

%% Reynolds - to be moved out of here

chord_eff = interp1(r(1,:),chord(1,:),0.7);
vel_eff = ((rotor.omega*rotor.R*0.7).^2+(interp1(r(1,:),lambda(1,:),0.7)).^2).^(1/2);
Re_eff = getReynolds(vel_eff,chord_eff,atm.kin_visc);

v_tip = rotor.omega*rotor.R;

velocity_dimensional = sqrt((v_tip*lambda_T.*sin(psi)+v_tip*r).^2+lambda.^2);

Re =  getReynolds(velocity_dimensional,chord,atm.kin_visc); % to be used in the future in get2Dcoeffs

%% Calculate dimensional parameters

Thrust = atm.rho*pi*rotor.R^4*rotor.omega^2*CT; %Using rotor 1 radius=rotor2 radius
Power = atm.rho*pi*rotor.R^5*rotor.omega^3*CP; %Using rotor 1 radius=rotor2 radius
Torque = atm.rho*pi*rotor.R^5*rotor.omega^2*CP; %Using rotor 1 radius=rotor2 radius

[alpha_negatives,AoA] = getAoA(rotor.pitch,phi_negative);


end
