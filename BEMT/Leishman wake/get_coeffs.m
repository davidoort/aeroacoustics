function [CT_u, CP_u, CT_l, CP_l] = get_coeffs(Fcf_u, lambda_u, Fcf_l, lambda_l, r, dr, rotor)

pitch_u = rotor(1).pitch;
pitch_l = rotor(2).pitch;

%% UPPER ROTOR
dcT_u = 4*Fcf_u.*lambda_u.^2.*r*dr;
phi_u = lambda_u./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

dcp_i_u = lambda_u.*dcT_u; %induced power/drag
Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_u = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(pitch_u-phi_u).*r.^3)+sum(D2*dr*(pitch_u-phi_u).^2.*r.^3)); %profile power/drag

CT_u = sum(dcT_u);
CP_u = sum(dcp_i_u)+cp_p_u; %the same as CQ_u

%% BOTTOM ROTOR

dcT_l = 4*Fcf_l.*lambda_l.^2.*r*dr;
phi_l = lambda_l./r; %rad - induced inflow angle. small angle approximation for tangent(phi) = phi

dcp_i_l = lambda_l.*dcT_l; %induced power/drag
Cd0 = rotor(1).aero.Cd0;
D1 = rotor(1).aero.D1;
D2 = rotor(1).aero.D2;
cp_p_l = 0.5*rotor(1).solidity*(sum(Cd0*dr*r.^3)+sum(D1*dr*(pitch_l-phi_l).*r.^3)+sum(D2*dr*(pitch_l-phi_l).^2.*r.^3)); %profile power/drag

CT_l = sum(dcT_l);
CP_l = sum(dcp_i_l)+cp_p_l; %the same as CQ_l



end