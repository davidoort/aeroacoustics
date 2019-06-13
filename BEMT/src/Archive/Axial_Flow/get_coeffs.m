function [C_P, CT_u, CP_u, CT_l, CP_l, spanwise_coeffs] = get_coeffs(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, coaxial, flowfield)
%{
AXIAL Flight
This function is a simple code implementation of the equations found in
literature (see below) for the thrust and torque/power coefficients (and 
their distributions along the span) of the top and bottom rotor.

Inputs:
    Fcf_u - (array) containing the converged (for a given collective pitch
    setting) Prandtl tip loss function at every blade span for upper rotor

    lambda_u - (array) with converged inflow distribution as a function of
    blade span for the upper rotor (excludes lambda_inf)

    Fcf_l - (array) containing the converged (for a given collective pitch
    setting) Prandtl tip loss function at every blade span for upper rotor

    lambda_l - (array) with converged inflow distribution as a function of
    blade span for the upper rotor (excludes lambda_inf)

    r - (array) of radial positions from dr to 1-dr

    dr - (scalar) spanwise (radial) discretization size

    rotor - (struct) containing geometrical properties for both rotors such 
    as pitch distribution, radius, rpm, blade number, etc 

    params - (struct) containing general parameters such as kappa and
    kappa_int

    flowfield - (struct) containing lambda_inf (array) for both rotors 
    (normalization is different for each rotor since tip speed may be different)

Outputs:
    C_P - (scalar) power coefficient for coaxial rotor system, from Leishman 
    equation (20), verification purposes

    CT_u - (scalar) upper rotor thrust coefficient. From Leishman equation (15)
    
    CP_u - (scalar) upper rotor power coefficient. From Leishman equation (16)

    CT_l - (scalar) lower rotor thrust coefficient. From Leishman equation (17)
 
    CP_l - (scalar) lower rotor power coefficient. From Leishman equation (18)

    spanwise_coeffs - (struct) containing the spanwise thrust and power
    coefficients, for plotting.

Other m-files required: none

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (13) & (14)

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions: none

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
March 2019; Last revision: 7-June-2019
%}

%------------- BEGIN CODE --------------

rotor = coaxial.rotor;

params = coaxial.params;

pitch_u = rotor(1).pitch;
pitch_l = rotor(2).pitch;

%% UPPER ROTOR
dcT_u = 4*Fcf_u.*lambda_u.*(lambda_u+flowfield(1).lambda_inf).*r*dr;
phi_u = (lambda_u+flowfield(1).lambda_inf)./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

alpha_u = pitch_u-phi_u;

alpha_0_u = rotor(1).aero.alpha_0;

%Different versions presented in the paper. 1 is wrong! 2&3 are the same in
%hover but I don't know which is correct for axial flight
%dcp_i_u_1 = 4*(lambda_u+flowfield(1).lambda_inf).^2.*lambda_u.*r*dr; %induced power/drag as in eq.(5b)
dcp_i_u_2 = lambda_u.*dcT_u; %induced power/drag - including F (lower value) and multiplied by lambda_u as in eq.(16)
dcp_i_u_3 = (lambda_u+flowfield(1).lambda_inf).*dcT_u; %induced power/drag - including F (lower value) and multiplied by lambda_u as in eq.(5a)

Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_u = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(alpha_u-alpha_0_u).*r.^3)+sum(D2*dr*(alpha_u-alpha_0_u).^2.*r.^3)); %profile power/drag

CT_u = sum(dcT_u);
CP_u = sum(dcp_i_u_3)+cp_p_u; %the same as CQ_u


spanwise_coeffs.dCp_u = dcp_i_u_3/dr + 0.5*rotor(1).solidity*(Cd0*r.^3+D1*(alpha_u-alpha_0_u).*r.^3+D2*(alpha_u-alpha_0_u).^2.*r.^3);
spanwise_coeffs.dCt_u = dcT_u/dr;

%% BOTTOM ROTOR

%dcT_l=dcT_l_1 !!
dcT_l_1 = 4*Fcf_l.*lambda_l.*(lambda_l+flowfield(2).lambda_inf).*r*dr;
%Alternate approach, from Leishman (there is actually no difference between the two expressions, I have checked)
dcT_l = 0.5*rotor(2).solidity*rotor(2).aero.cl_alpha*(rotor(2).pitch.*r.^2-...
    (lambda_l+flowfield(2).lambda_inf).*r)*dr; %the SAME as what is written above

phi_l = (lambda_l+flowfield(2).lambda_inf)./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

alpha_l = pitch_l-phi_l;
alpha_0_l = rotor(2).aero.alpha_0;
%Different versions presented in the paper. 1 is wrong! 2&3 are the same in
%hover but I don't know which is correct for axial flight
%dcp_i_l_1 = 4*(lambda_l+flowfield(2).lambda_inf).^2.*lambda_l.*r*dr; %induced power/drag as in eq.(5b)
dcp_i_l_2 = lambda_l.*dcT_l; %induced power/drag - including F (lower value) and multiplied by lambda_u as in eq.(16)
dcp_i_l_3 = (lambda_l+flowfield(2).lambda_inf).*dcT_l; %induced power/drag - including F (lower value) and multiplied by lambda_u as in eq.(5a)

Cd0 = rotor(2).aero.Cd0;
D1 = rotor(2).aero.D1;
D2 = rotor(2).aero.D2;
cp_p_l = 0.5*rotor(2).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(alpha_l-alpha_0_l).*r.^3)+sum(D2*dr*(alpha_l-alpha_0_l).^2.*r.^3)); %profile power/drag

CT_l = sum(dcT_l);
CP_l = sum(dcp_i_l_3)+cp_p_l; %the same as CQ_l


spanwise_coeffs.dCp_l = dcp_i_l_3/dr + 0.5*rotor(2).solidity*(Cd0*r.^3+D1*(alpha_l-alpha_0_l).*r.^3+D2*(alpha_l-alpha_0_l).^2.*r.^3);
spanwise_coeffs.dCt_l = dcT_l/dr;


%% C_P New
%In the Leishman paper it said that using the BEMT he could calculate CP
%and CT only by estimating Cd0, so the above formula has to be replaced by:

C_P = CP_u+CP_l;


end