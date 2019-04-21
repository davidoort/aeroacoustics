function [output1, output2] = templatefunc(input1, input2, input3)
%{
templatefunc - General description

Inputs:
    input1 - (data type) Description

    input2 - (data type) Description
    
    input3 - (data type) Description

Outputs:
    output1 - (data type) Description

    output2 - (data type) Description

Other m-files (functions) required (added to path): 
    get_lambda_up
    Prandtl_tip_loss
    get_lambda_bot

MAT-files required: none

Literature used: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008.

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions:
    ! Top rotor not influenced by bottom rotor: its inflow can be converged
    independently.

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
March 2019; Last revision: 21-April-2019
%}

%------------- BEGIN CODE --------------

disp("template")

%------------- END OF CODE --------------

end