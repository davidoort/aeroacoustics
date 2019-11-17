function lambda = getLambda_Leish(rotor,lambda_ext,lambda_T,F,r,psi)
%{
getLambda_Leish APPROACH1 finds the induced velocity of a (or multiple) disk element(s)
when given an external axial vel (which could also include the downwash of 
another rotor), a prandtl tip loss function

Inputs:
    lambda_ext - (matrix [-]) lambda_P for upper rotor and
    lambda_P+wake_upper for lower rotor affected area

    F - (matrix [-]) Prandtl tip loss function for each element

    r - (matrix) non-dimensional radial position to be examined

    dr - (scalar) width of non-dimensional radial discretization

    dpsi - (scalar [rad]) width of dimensional azimuthal discretization

Outputs:

    dCT - (scalar, array or matrix [-]) incremental thrust coefficient of annulus/annuli
    
Other m-files required: none

Other files required: none

Literature referenced: 

    An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. 

    ETH slides - for basics of momentum theory

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 13-June-2019
%}

%------------- BEGIN CODE --------------

mu = lambda_T.^2.*(sin(psi)).^2+2*lambda_T.*r.*sin(psi);

b = rotor.solidity*rotor.aero.cl_alpha*(mu+r.^2)./(16*F.*r.*(r+lambda_T.*sin(psi)))-lambda_ext/2;

lambda = sqrt(b.^2+rotor.solidity*rotor.aero.cl_alpha*rotor.pitch.*(mu+r.^2)./(8*F.*r))-b;

end