function [CT_u, CP_u, spanwise_coeffs] = get_coeffs_single(Fcf_u, lambda_u, r, dr, single, flowfield)
%{
This function is a simple code implementation of the equations found in
literature (see below) for the thrust and torque/power coefficients (and 
their distributions along the span) of the top (SINGLE) rotor.

Inputs:
    Fcf_u - (array) containing the converged (for a given collective pitch
    setting) Prandtl tip loss function at every blade span for upper rotor

    lambda_u - (array) with converged inflow distribution as a function of
    blade span for the upper rotor

    r - (array) of radial positions from dr to 1-dr

    dr - (scalar) spanwise (radial) discretization size

    rotor - (struct) containing geometrical properties for both rotors such 
    as pitch distribution, radius, rpm, blade number, etc 

    params - (struct) containing general parameters such as kappa and
    kappa_int

    flowfield - (struct) containing lambda_inf (array) for both rotors 
    (normalization is different for each rotor since tip speed may be different)

Outputs:

    CT_u - (scalar) upper rotor thrust coefficient. From Leishman equation (15)
    
    CP_u - (scalar) upper rotor power coefficient. From Leishman equation (16)

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
May 2019; Last revision: 8-May-2019
%}

%------------- BEGIN CODE --------------

rotor = single.rotor(1);

pitch_u = rotor(1).pitch;

%% UPPER ROTOR
dcT_u = 4*Fcf_u.*lambda_u.*(lambda_u+flowfield(1).lambda_inf).*r*dr;
phi_u = (lambda_u+flowfield(1).lambda_inf)./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

alpha_u = pitch_u-phi_u;
alpha_0_u = rotor(1).aero.alpha_0;

dcp_i_u = lambda_u.*dcT_u; %induced power/drag
Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_u = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(alpha_u-alpha_0_u).*r.^3)+sum(D2*dr*(alpha_u-alpha_0_u).^2.*r.^3)); %profile power/drag

CT_u = sum(dcT_u);
CP_u = sum(dcp_i_u)+cp_p_u; %the same as CQ_u


spanwise_coeffs.dCp_u = dcp_i_u/dr + 0.5*rotor(1).solidity*(Cd0*r.^3+D1*(alpha_u-alpha_0_u).*r.^3+D2*(alpha_u-alpha_0_u).^2.*r.^3);
spanwise_coeffs.dCt_u = dcT_u/dr;


end