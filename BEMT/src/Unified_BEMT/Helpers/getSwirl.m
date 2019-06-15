function [weighted_swirl_ratio] = getSwirl(lambda,lambda_ext,lambda_T,r,dr,dpsi,dCP)
%{
getSwirl is a helper function that calculates the rotational speed of the
wake of a rotor normalized by its angular velocity.

PROBABLY DOESN'T WORK IN FORWARD FLIGHT - should be extended at least

Inputs:

    r - (array, scalar or matrix) non-dimensional radial position to be examined

    dr - (scalar) width of non-dimensional radial discretization

    lambda - (array or scalar (matrix?)) non-dimensional inflow ratio (perpendicular to blade) at each radial position
    
    dCP - (array) power coefficient of the rotor at each annulus. A bit
    confused about whether this needs to be multiplied by dr or not. My
    guess is that it should, for the equation with a dr in the denominator.

Outputs:

    weighted_swirl_ratio - (scalar [-]) a massflow-weighted swirl ratio of
    the propeller wake. 
    
Other m-files required: none

Other files required: none

Literature referenced: 
    
    BEM theory and CFD for Wind Turbine Aerodynamics. Internship report

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------


%swirl ratio: swirl wake/omega

%mass_flow_element/(rho*V_tip) = dA*lambda = lambda*r*dr*dpsi

if norm(lambda_T)==0 %otherwise it might not be a reasonable assumption
    lambda(isnan(lambda))=lambda_ext(1,1); 
end

mass_flow_element = lambda.*r*dr*dpsi;

if sum(sum(mass_flow_element)) == 0
    disp('No mass flow through the rotor')
    weighted_swirl_ratio = 0;
else
    swirl_ratio = pi*dCP./(mass_flow_element.*r.^2);

    weighted_swirl_ratio = sum(sum(swirl_ratio.*mass_flow_element))/sum(sum(mass_flow_element));
end

end

