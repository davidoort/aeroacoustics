function lambda = getLambda(lambda_ext,dCT,F,r,dr,dpsi)
%{
getLambda finds the induced velocity of a (or multiple) disk element(s)
when given an external axial vel (which could also include the downwash of 
another rotor), a thrust coefficient and a prandtl tip loss function

Inputs:
    lambda_ext - (matrix [-]) lambda_P for upper rotor and
    lambda_P+wake_upper for lower rotor affected area

    dCT - (matrix [-]) thrust coefficient of each element

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
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------

negCT_bool = dCT<0;
posCT_bool = dCT >= 0;

lambda_i_pos = sqrt((lambda_ext/2).^2+pi*dCT./(2*F.*r*dr*dpsi))-(lambda_ext/2); %normal formula



lambda_up = sqrt((lambda_ext/2).^2+pi*-dCT./(2*F.*r*dr*dpsi))+(lambda_ext/2); %lambda up as defined in my notes



lambda_i = lambda_i_pos.*posCT_bool-lambda_up.*negCT_bool;



lambda = lambda_i +lambda_ext;

end

