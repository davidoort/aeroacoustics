function [Cl,Cd,Cl_max,alpha_stall,Cds] = lookup_coeffs(Re,alpha,airfoil)
%{
LOOKUP_COEFFS is a helper function that returns aerodynamic coefficients
for a specified airfoil, Reynolds number and angle of attack based on a lookup table from XFOIL. 

Limited by the fact that the only look up tables are for NACA 0012 and do
not include boundary layer properties.

Inputs:
    Re - (scalar) Reynolds number (note that instead I could have given it 
    atm and other rotor properties but the coding philosophy here is to let
    this function be very light and transparent)

    alpha - (scalar [deg]) with atmospheric parameters such as air density

    airfoil - (string) name of the airfoil

Outputs:

    Cl - (matrix) non-dimensional 2D lift coefficient of a blade element
    
    Cd - (matrix) non-dimensional 2D drag coefficient of a blade element

Other m-files required: none

Other files required: 

    Airfoil database

Literature referenced: 

    XFOIL

    Airfoiltools - http://airfoiltools.com/airfoil/details?airfoil=n0012-il

Assumptions: Linear interpolation
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 10-June-2019
%}

%------------- BEGIN CODE --------------

Re = num2str(Re);
airfoil = string(airfoil);

%% Sanity checks


if strcmpi(airfoil,'NACA0012')
    data = readmatrix('xf-n0012-il-1000000'); % Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr
elseif strcmpi(airfoil,'NACA16006')
    data = readmatrix('xf-naca16006-il-1000000'); % Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr
elseif strcmpi(airfoil,'ClarkY')
    data = readmatrix('xf-clarky-il-1000000-n5'); % Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr
elseif strcmpi(airfoil,'NACA23015')
    data = readmatrix('xf-naca23015-il-1000000'); % Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr
else
    error('Currently there is only data for NACA0012, NACA16006 and ClarkY (Re=1e6). Please download the airfoil data of whatever airfoil you are using')
end


%ideally, this only shows up with verbose or debug
% if max(alpha)>18.5
%     warning(strcat('You have passed an alpha of ', num2str(max(alpha)), ' > 18.5 deg. Using high angle of attack corrections...'))
% elseif min(alpha)<-18.5
%     warning(strcat('You have passed an alpha of ', num2str(min(alpha)), ' < -18.5 deg. Using high angle of attack corrections...'))
% end

alpha_arr = data(:,1);
Cl_arr = data(:,2);
Cd_arr = data(:,3);

%Cl_v_alpha = griddedInterpolant(alpha_arr,Cl_arr);
%Cl = Cl_v_alpha(alpha);

%griddedInterpolant might be faster

Cl_max = max(Cl_arr);

alpha_stall = interp1(Cl_arr,alpha_arr,Cl_max);

Cl = interp1(alpha_arr,Cl_arr,alpha);
Cd = interp1(alpha_arr,Cd_arr,alpha);
Cds = interp1(alpha_arr,Cd_arr,alpha_stall);


end

