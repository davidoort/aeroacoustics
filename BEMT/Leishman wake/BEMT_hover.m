
%maybe make this a function, like get_coeffs(collective1,collective2). This
%allows to trim the helicopter in hover. Whenever you change pitch and
%therefore thrust and torque it would be good to start the convergence of
%Fcf and lambda (in a separate function) with the values previously found.

init 
dr = 0.001;
r = dr:dr:1-dr; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

epsilon = 0.001; %convergence accuracy

pitch_l = 2/12*ones(1,length(r));
pitch_u = pitch_l; %assumed that all blades (on both rotors) have the same geometric pitch and are at the same collective setting

%% Converge Fcf and inflow ratio 

    
%% Calculate thrust and torque coefficients

%UPPER ROTOR
dcT_u = 4*Fcf_u.*lambda_u.^2.*r*dr;
phi_u = lambda_u./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

dcp_i_u = lambda_u.*dcT_u; %induced power/drag
Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_u = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(pitch_u-phi_u).*r.^3)+sum(D2*dr*(pitch_u-phi_u).^2.*r.^3)); %profile power/drag

CT_u = sum(dcT_u);
CP_u = sum(dcp_i_u)+cp_p_u; %the same as CQ_u
FOM_u = CT_u^(3/2)/(sqrt(2)*CP_u); %treated as a single rotor;

%BOTTOM ROTOR

dcT_l = 4*Fcf_l.*lambda_l.^2.*r*dr;
phi_l = lambda_l./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

dcp_i_l = lambda_l.*dcT_l; %induced power/drag
Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_l = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(pitch_l-phi_l).*r.^3)+sum(D2*dr*(pitch_l-phi_l).^2.*r.^3)); %profile power/drag

CT_l = sum(dcT_l);
CP_l = sum(dcp_i_l)+cp_p_l; %the same as CQ_l
FOM_l = CT_l^(3/2)/(sqrt(2)*CP_l); %treated as a single rotor;


