function [Fcf] = getPrandtlTipLoss(rotor,phi,r)
%{
getPrandtlTipLoss is a code implementation of the equation found in
literature (see below) for the Prandtl Tip Loss function F as a function of
blade number and spanwise inflow ratio. Function performs calculation in a
disk-element-wise fashion.

Inputs:
    lambda - (array) of inflow distribution (not necessarily converged)
    along the span.

    r - (array) 

    rotor - (struct) containing geometrical properties for ONE of the two rotors such 
    as pitch distribution, radius, rpm, blade number, etc. Should be passed
    either rotor(1) or rotor(2)

Outputs:
    Fcf - (array) Prandtl tip loss function (not necessarily converged) 
    as a function of blade number and inflow.

Other m-files required: none

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (6) & (7)

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions:
    REMOVED SMALL ANGLE ASSUMPTION

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}



%------------- BEGIN CODE --------------


Fcf_tip = (2/pi)*acos(exp((-rotor.Nb/2)*(ones(size(r))-r)./(r.*phi)));
Fcf_root = (2/pi)*acos(exp((-rotor.Nb/2)*(r-rotor.hub_radial_fraction*ones(size(r)))./(r.*phi)));


Fcf = Fcf_tip.*Fcf_root; %doesn't make much of a difference in the inflow plot, maybe I need a different F function (like the one from Carlos..)



end

