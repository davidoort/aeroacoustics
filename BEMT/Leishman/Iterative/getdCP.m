function dCP = getdCP(rotor,atm,phi,r,dr,psi,dpsi,chord,lambda,lambda_T)
%{
GETdCP calculates the thrust coefficient of one
(or multiple) rotor disk elements

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

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------

[Cl,Cd] = get2Dcoeffs(rotor,atm,phi,chord,lambda,lambda_T,r,psi);
                       
dCP = (rotor.Nb*chord.*r*dr)/(2*pi*rotor.R).*((lambda_T.*sin(psi)+r).^2+lambda.^2).*(Cl.*sin(phi)+Cd.*cos(phi))*dpsi/(2*pi);

end

