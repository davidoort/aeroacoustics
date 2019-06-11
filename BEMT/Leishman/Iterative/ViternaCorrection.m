function [Cl,Cd] = ViternaCorrection(aspect_ratio, Cl_orig, Cd_orig, alpha, alpha_s, Cds, Clmax)
%{
ViternaCorrection is a helper function that calculates lift and drag
coefficients in stall

Inputs:

    aspect_ratio - (scalar) aspect ratio of blade

    Cl_orig - (matrix [-]) original matrix of Cl values at disk elements with
    nan entries for the angles of attack that need to be extrapolated from
    the lookup tables

    Cd_orig - (matrix [-]) original matrix of Cd values at disk elements with
    nan entries for the angles of attack that need to be extrapolated from
    the lookup tables

    alpha - (matrix [deg]) matrix of angles of attack at each disk element

    alpha_s - (scalar [deg]) stall angle of attack of the rotor airfoil

    Cds - (scalar [-]) drag coefficient from lookup table evaluated at alpha_s

    Clmax - (scalar [-]) maximum lift coefficient of the rotor airfoil (or Cls)

Outputs:

    Cl - (matrix) non-dimensional 2D lift coefficient of a blade element
    with nan entries corrected according to viterna

    Cd - (matrix) non-dimensional 2D drag coefficient of a blade element
    with nan entries corrected according to viterna
    
Other m-files required: none

Other files required: none

Literature referenced: 

    Fixed Pitch Rotor Performance of Large Horizontal Axis Wind Turbines.
    Larry A. Viterna & Robert D. Corrigan.  NASA CP 2230, DOE Publication 
    CONF 810752, Cleveland, OH, NASA Lewis Research Centre, 1981.
    
    Improving BEM‐based Aerodynamic Models in Wind Turbine Design Codes. PhD
    Thesis Tonio Sant

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------


Cl_stalled_bool = isnan(Cl_orig);
Cd_stalled_bool = isnan(Cd_orig);

Cl_unstalled_bool = ~isnan(Cl_orig);
Cd_unstalled_bool = ~isnan(Cd_orig);

if aspect_ratio<=50 && aspect_ratio>0
    Cdmax = 1.11+0.018*aspect_ratio;
else
    Cdmax = 2.01;
end

K_l = (Clmax-Cdmax*sin(deg2rad(alpha_s))*cos(deg2rad(alpha_s)))*sin(deg2rad(alpha_s))/(cos(deg2rad(alpha_s)))^2;
K_d = (Cds-Cdmax*(sin(deg2rad(alpha_s)))^2)/cos(deg2rad(alpha_s));

Cl_viterna = Cdmax*sin(2*deg2rad(alpha))/2+K_l*cos(deg2rad(alpha))^2/sin(deg2rad(alpha));
Cd_viterna = Cdmax*(sin(deg2rad(alpha)))^2+K_d*cos(deg2rad(alpha));

Cl = Cl_viterna.*Cl_stalled_bool + Cl_orig.*Cl_unstalled_bool;
Cd = Cd_viterna.*Cd_stalled_bool + Cd_orig.*Cd_unstalled_bool;


end

