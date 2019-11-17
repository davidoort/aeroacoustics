function [Cl,Cd] = get2Dcoeffs(rotor,atm,phi,chord,lambda,lambda_T,r,psi)
%{
get2Dcoeffs is a helper function that returns aerodynamic coefficients
for a specified airfoil, Reynolds number and angle of attack based on a
lookup table from XFOIL and Viterna high angle of attack corrections

Inputs:
    
    rotor - (struct) geometric variables (such as airfoil) of A rotor

    atm - (struct) struct containing atmospheric parameters such as air
    density, kinematic viscosity, etc.

    phi - (matrix [rad]) inflow angle

    chord - (matrix [m]) chord of the rotor blade at each radial and
    azimuthal location

    lambda - (matrix) non-dimensional inflow ratio (perpendicular to blade) at each radial position
    
    lambda_T - (scalar) non-dimensional inflow ratio (tangent to blade) at
    each radial position 
    
    psi - (matrix [rad]) dimensional azimuthal position to be examined

    r - (matrix) non-dimensional radial position to be examined

Outputs:

    Cl - (matrix) non-dimensional 2D lift coefficient of a blade element
    
    Cd - (matrix) non-dimensional 2D drag coefficient of a blade element

Other m-files required: none

Other files required: 

    getReynolds

    lookup_coeffs

    getAspectRatio

    ViternaCorrection

Literature referenced: In the functions called

Assumptions: Linear interpolation
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 11-June-2019
%}

%------------- BEGIN CODE --------------

%chord_eff = interp1(r(1,:),chord(1,:),0.7);
%vel_eff = ((rotor.omega*rotor.R*0.7).^2+(interp1(r(1,:),lambda(1,:),0.7)).^2).^(1/2);
%Re_eff = getReynolds(vel_eff,chord_eff,atm.kin_visc);

v_tip = rotor.omega*rotor.R;

velocity_dimensional = sqrt((v_tip*lambda_T.*sin(psi)+v_tip*r).^2+lambda.^2);

Re =  getReynolds(velocity_dimensional,chord,atm.kin_visc); % for later, when more lookup tables exist

airfoil = rotor.airfoil.name;

[AoA,AoA_nans] = getAoA(rotor.pitch,phi);
[Cl,Cd,Cl_max,alpha_stall,Cds] = lookup_coeffs(1e6,AoA_nans,airfoil); %the 1e6 is a quick bug fix

negthrust_bool = isnan(AoA_nans); % here assuming that positive alpha is also positive Cl which is true for normal airfoils that don't have negative lift at alpha 0

Cl_positive = Cl; 
Cl_positive(negthrust_bool) = 0; %to be passed to Viterna 
Cd_for_Cl_positive = Cd; 
Cd_for_Cl_positive(negthrust_bool) = 0; %to be passed to Viterna

stalled_bool = isnan(Cl_positive);
not_stalled_bool = ~isnan(Cl_positive);

if any(any(isnan(Cl_positive))) || any(any(isnan(Cd_for_Cl_positive))) %Cl_positive should have only nans remaining for high AoA
    hub_R = rotor.hub_radial_fraction*rotor.R;
    aspect_ratio = getAspectRatio(rotor.root_chord,rotor.tip_chord,hub_R,rotor.R);
    [Cl_positive,Cd_for_Cl_positive] = ViternaCorrection(aspect_ratio, Cl_positive, Cd_for_Cl_positive, AoA, alpha_stall, Cds, Cl_max);    
end

%Cl(~isnan(AoA_nans)) = Cl_positive; %reshaped array, for some reason, better to use logical matrices
%Cd(~isnan(AoA_nans)) = Cd_for_Cl_positive;
Cl(stalled_bool) = 0;
Cd(stalled_bool) = 0;

Cl = Cl.*not_stalled_bool+Cl_positive.*stalled_bool;
Cd = Cd.*not_stalled_bool+Cd_for_Cl_positive.*stalled_bool;


end

