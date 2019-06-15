function dCT = getdCT_Leish(rotor,F,lambda,lambda_P,lambda_T,phi,r,dr,psi,dpsi)
%{
GETdCT_Leish calculates the thrust coefficient of one
(or multiple) rotor disk elements as derived by me in eq (21)

Inputs:
    rotor - (struct) geometric variables (such as airfoil) of A rotor

    atm - (struct) struct containing atmospheric parameters such as air
    density, kinematic viscosity, etc.

    phi - (matrix [rad]) inflow angle

    r - (matrix) non-dimensional radial position to be examined

    dr - (scalar) width of non-dimensional radial discretization

    lambda - (matrix) non-dimensional inflow ratio (perpendicular to blade) at each radial position
    
    lambda_T - (scalar) non-dimensional inflow ratio (tangent to blade) at
    each radial position 
    
    psi - (matrix [rad]) dimensional azimuthal position to be examined

    dpsi - (scalar [rad]) width of dimensional azimuthal discretization

Outputs:

    dCT - (scalar, array or matrix [-]) incremental thrust coefficient of annulus/annuli
    
Other m-files required: none

Other files required: none

Literature referenced: 

    ETH slides - basics of momentum theory

    Honours report

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 13-June-2019
%}

%------------- BEGIN CODE --------------


alpha = rotor.pitch-phi;

alpha_0 = rotor(1).aero.alpha_0;

dCT1 = rotor.solidity*rotor.aero.cl_alpha*(alpha-alpha_0).*(lambda_T.^2.*sin(psi).^2+2*lambda_T.*r.*sin(psi)+r.^2)*dr*dpsi/(4*pi); % equation 16
dCT = 2*F.*lambda.*(lambda-lambda_P).*r*dr*dpsi/pi; %equation 23

end

