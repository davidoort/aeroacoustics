function lambda_bot = get_lambda_bot(F,r,pitch,rotor)

%F = Prandtl tip loss
%r = radial position of blade
%pitch of blade element wrt to rotor plane in rad
%rotor = rotor struct

rotor_full = rotor; %used to call lambda_up
rotor = rotor(2); %in case the generic rotor struct is given
sigma = rotor.solidity;
cl_a = rotor.cl_alpha;

if r<rotor.rd
    lambda_up = get_lambda_up(F,r,pitch,rotor_full);
    lambda_bot = sqrt((sigma*cl_a/(16*F)-lambda_up/(2*rotor.rd^2))^2+sigma*cl_a*pitch*r/(8*F))-sigma*cl_a/(16*F)+lambda_up/(2*rotor.rd^2);
else
    lambda_bot = sqrt((sigma*cl_a/(16*F))^2+sigma*cl_a*pitch*r/(8*F))-sigma*cl_a/(16*F);
end







end