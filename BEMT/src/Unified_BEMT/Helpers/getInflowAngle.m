function [phi,phi_nans] = getInflowAngle(lambda,lambda_P,r,psi,lambda_T)
%{
GETINFLOWANGLE is a helper function that returns the inflow angle at a radial
station of the rotor. 

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

tangential_vel = (r+lambda_T.*sin(psi));

%tangential_vel = r;

backflow_bool = tangential_vel<=0;

lambda_0s = lambda;
lambda_0s(isnan(lambda))=lambda_P(isnan(lambda));

phi = atan(lambda_0s./tangential_vel);

phi_nans = atan(lambda./tangential_vel);
phi_nans(backflow_bool)=nan; %backflow

end

