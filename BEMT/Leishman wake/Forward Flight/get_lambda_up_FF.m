function lambda = get_lambda_up_FF(F,r,phi,pitch,rotor,flowfield)
%{
FORWARD Flight
This function is a simple code implementation of the equation found in
literature (see below) for the inflow distribution of the top rotor.

Inputs:
    F - (array) Prandtl tip loss function (not necessarily converged) 
    on the top rotor. 

    r - (array) of radial positions from dr to 1-dr

    pitch - (array) of pitch of each of the blade elements with respect to
    the rotor plane. This will be influenced by geometrical twist of the
    blade as well as by the collective setting.

    rotor - (struct) containing geometrical properties for both rotors such 
    as pitch distribution, radius, rpm, blade number, etc 

    flowfield - (struct) containing lambda_P (array) for both rotors 
    (normalization is different for each rotor since tip speed may be different)

Outputs:
    lambda - (array) of inflow distribution (not necessarily converged)
    along the span of the top rotor.

Other m-files required: none

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (12)

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
May 2019; Last revision: 26-May-2019
%}

%------------- BEGIN CODE --------------


rotor = rotor(1); %in case the generic rotor struct is given

lambda_P = flowfield(1).lambda_P;
sigma = rotor.solidity;
cl_a = rotor.aero.cl_alpha;

lambda = sqrt((sigma*cl_a*1./(16*F)-lambda_P/2).^2+sigma*cl_a*pitch.*r*1./(8*F))-sigma*cl_a*1./(16*F)+lambda_P/2;



end