function [dCP] = getdCP_Leish(rotor,dCT,r,dr,psi,dpsi,lambda,lambda_P,lambda_T,Cd)
%{
GETdCP_Leish calculates the thrust coefficient of one
(or multiple) rotor disk elements using my own derivation from Leishman
method

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

lambda_i = lambda-lambda_P;

dCP_i = lambda.*dCT; %if I am correct about the assumption that it should be lambda and not lambda_i!!!! - testing both
%dCP_i = 0;

mu = lambda_T.^2.*(sin(psi)).^2 + 2*lambda_T.*r.*sin(psi); %probably should make this a function at some point

dCP_p = rotor.solidity*dr*dpsi*r.*Cd.*(mu+r.^2)/(4*pi);

dCP =dCP_i+dCP_p;


end

