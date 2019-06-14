function [AoA,AoA_nans] = getAoA(pitch,phi)
%{
getAoA is a helper function that calculates the angle of attack at a disk
element 

Inputs:
    lambda - (matrix [-]) non-dimensional inflow ratio

    r - (matrix [-]) non-dimensional radial location of the root
    
    psi - (matrix [rad]) matrix of azimuthal locations

    lambda_T - (scalar [-]) non-dimensionalized tangential velocity
Outputs:

    phi - (array or scalar [rad]) dimensional inflow angle at r
    
Other m-files required: none

Other files required: none

Literature referenced: none

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 10-June-2019
%}

%------------- BEGIN CODE --------------

AoA = rad2deg(pitch-phi); %kind of works like this, when cyclic is used then pitch will be a matrix


AoA_nans = AoA;
AoA_nans(AoA<0)=nan;

end

