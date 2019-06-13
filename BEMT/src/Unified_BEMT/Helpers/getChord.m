function [chord] = getChord(r,r_root,tip_chord,root_chord)
%{
GETCHORD is a helper function that returns the chord at a spanwise station of the blade 

Inputs:
    r - (array or scalar [-]) non-dimensional radial query location

    r_root - (scalar [-]) non-dimensional radial location of the root

    tip_chord - (scalar [m]) dimensional chord length at the tip

    root_chord - (scalar [m]) dimensional chord length at the root

Outputs:

    chord - (array or scalar [m]) dimensional chord length at r
    
Other m-files required: none

Other files required: none

Literature referenced: none

Assumptions: Linear taper
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 10-June-2019
%}

%------------- BEGIN CODE --------------

query_num = length(r);

chord = ones(1,query_num)*root_chord-(root_chord-tip_chord)/(1-r_root) .* (r-ones(1,query_num)*r_root); %m


end

