
%maybe make this a function, like get_coeffs(collective1,collective2). This
%allows to trim the helicopter in hover. Whenever you change pitch and
%therefore thrust and torque it would be good to start the convergence of
%Fcf and lambda (in a separate function) with the values previously found.

init 
dr = 0.01;
r = dr:dr:1-dr; %non-dimensionalized by tip radius. Rotors have the same radius.
%The dr and 1-dr is to avoid singularities at the tip (since F= 0 there usually and the lambda is NaN)
%and at the root, when calculating the induced inflow angle.

epsilon = 0.001; %convergence accuracy

pitch_l = 1/12*ones(1,length(r));
pitch_u = pitch_l; %assumed that all blades (on both rotors) have the same geometric pitch and are at the same collective setting

%% Converge Fcf and inflow ratio 

Fcf0_u = zeros(1,length(r)); %dummy, to start iteration
Fcf0_l = zeros(1,length(r)); %dummy, to start iteration
Fcf_u = ones(1,length(r)); %"Fcf started with an initial value of 1"
Fcf_l = ones(1,length(r)); %"Fcf started with an initial value of 1"

lambda0_u = zeros(1,length(r)); %dummy, to start iteration
lambda0_l = zeros(1,length(r)); %dummy, to start iteration
lambda_u = get_lambda_up(Fcf_u,r,pitch_u,rotor);
lambda_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor);


%converge upper rotor - can be done without converging bottom rotor (assumption of no influence of bottom on top)
i = 0;
while norm(Fcf_u-Fcf0_u)>epsilon || norm(lambda_u-lambda0_u)>epsilon
    Fcf0_u = Fcf_u;
    lambda0_u = lambda_u;
    
    Fcf_u = Prandtl_tip_loss(r,lambda_u,rotor(1));
    lambda_u = get_lambda_up(Fcf_u,r,pitch_u,rotor);
    i = i+1
end

%converge bottom rotor 
i = 0;
while norm(Fcf_l-Fcf0_l)>epsilon || norm(lambda_l-lambda0_l)>epsilon
    Fcf0_l = Fcf_l;
    lambda0_l = lambda_l;
    
    Fcf_l = Prandtl_tip_loss(r,lambda_l,rotor(2));
    lambda_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor);
    i = i+1
end
    
    
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
CP_u = sum(dcp_i_u)+cp_p_u;

%BOTTOM ROTOR


%SPECS

%FOM = 

