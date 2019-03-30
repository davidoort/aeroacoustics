function lambda_up = get_lambda_up(F,r,pitch,rotor)

%F = Prandtl tip loss
%r = radial position of aero
%pitch of blade element wrt to rotor plane in rad
%rotor= rotor struct

rotor = rotor(1); %in case the generic rotor struct is given

sigma = rotor.solidity;
cl_a = rotor.aero.cl_alpha;

lambda_up = sqrt((sigma*cl_a*1./(16*F)).^2+sigma*cl_a*pitch.*r*1./(8*F))-sigma*cl_a*1./(16*F);



end