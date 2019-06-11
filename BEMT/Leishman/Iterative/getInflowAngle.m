function [phi] = getInflowAngle(lambda,r)
%{
GETINFLOWANGLE is a helper function that returns the inflow angle at a radial
station of the rotor. 
MIGHT GET EXTENDED TO FORWARD FLIGHT IN THE SPIRIT OF UNITING THE CODE

Inputs:
    lambda - (array or scalar [-]) non-dimensional inflow ratio

    r - (array or scalar [-]) non-dimensional radial location of the root

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

phi = atan(lambda./r);

end

