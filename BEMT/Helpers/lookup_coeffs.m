function [Cl,Cd] = lookup_coeffs(Re,alpha,airfoil)
%{
LOOKUP_COEFFS is a helper function that returns aerodynamic coefficients
for a specified airfoil, Reynolds number and angle of attack. 

Limited by the fact that the only look up tables are for NACA 0012 and do
not include boundary layer properties.

Inputs:
    Re - (scalar) Reynolds number (note that instead I could have given it 
    atm and other rotor properties but the coding philosophy here is to let
    this function be very light and transparent)

    alpha - (struct object) with atmospheric parameters such as air density

    airfoil - (string) name of the airfoil

    CT_or_pitch - (scalar) can be either collective_u (upper rotor pitch 
    angle in degrees) or the desired thrust coefficient of the coaxial
    system

    trimvar - (string) indicates if CT desired or pitchdeg desired is being 
    specified

Outputs:

    Cl - (scalar) non-dimensional 2D lift coefficient of a blade element
    
    Cd - (scalar) non-dimensional 2D drag coefficient of a blade element

Other m-files required: none

Other files required: 

    Airfoil database

Literature referenced: 

    XFOIL

    Airfoiltools - http://airfoiltools.com/airfoil/details?airfoil=n0012-il

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 8-June-2019
%}

%------------- BEGIN CODE --------------

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

