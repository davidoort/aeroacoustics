function [Re] = getReynolds(vel,chord,kin_visc)
%{
getReynolds is a helper function that returns the chord at a spanwise station of the blade 

Inputs:
    
    vel - (scalar or array [m/s]) velocity of the fluid wrt to body

    chord - (scalar or array [m]) characteristic linear dimension which for blades
    is the chord

    kin_visc - (scalar [m^2/s]) kinematic viscosity of the fluid

Outputs:

    Re - (scalar or array [-]) Reynolds number
    
Other m-files required: none

Other files required: none

Literature referenced: Airfoiltools, Wikipedia

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 10-June-2019
%}

%------------- BEGIN CODE --------------

Re = vel.*chord/kin_visc;

end

