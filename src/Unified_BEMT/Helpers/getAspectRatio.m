function aspect_ratio = getAspectRatio(root_chord,tip_chord,hub_R,R)
%{
getAspectRatio is a helper function that calculates the aspect ratio of a
LINEARLY TAPERED blade

Inputs:

    root_chord - (scalar [m]) dimensional chord length at the root

    tip_chord - (scalar [m]) dimensional chord length at the tip

    hub_R - (scalar [m]) Radius of the hub

    R - (scalar [m]) Radius of the rotor blade

Outputs:

    aspect_ratio - (scalar [-]) aspect ratio of blade (same definition as that for a wing)
    
Other m-files required: none

Other files required: none

Literature referenced: 

    Rule of thumb b^2/S

Assumptions: LINEAR TAPER
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------

blade_length = R-hub_R;

S = tip_chord*blade_length+(tip_chord-root_chord)*blade_length/2;

aspect_ratio = blade_length^2/S;

end

