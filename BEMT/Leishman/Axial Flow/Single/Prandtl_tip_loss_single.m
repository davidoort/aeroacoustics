function Fcf= Prandtl_tip_loss_single(r,lambda,rotor)
%{
This function is a simple code implementation of the equation found in
literature (see below) for the Prandtl Tip Loss function F as a function of
blade number and spanwise inflow ratio. Function performs calculation in a
blade-element-wise fashion.

Inputs:
    lambda - (array) of inflow distribution (not necessarily converged)
    along the span.

    r - (array) of radial positions from dr to 1-dr

    rotor - (struct) containing geometrical properties for ONE of the two rotors such 
    as pitch distribution, radius, rpm, blade number, etc. Should be passed
    either rotor(1)

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
    ! Small angle assumption being used in equation 7

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
May 2019; Last revision: 8-May-2019
%}



%------------- BEGIN CODE --------------

dim = size(rotor);

if dim(2)>1
    error("Function expects either rotor(1) or rotor(2)")
end

Fcf = (2/pi)*acos(exp(-rotor.blades.num*(ones(1,length(r))-r)./(2*lambda)));

%------------- END OF CODE --------------

end
