function lambda_bot = get_lambda_bot(F,r,pitch,rotor,flowfield)

%F = Prandtl tip loss (array)
%r = radial position of blade (array)
%pitch of blade element wrt to rotor plane in rad (array)
%rotor = rotor struct

lambda_inf = flowfield(2).lambda_inf;

rotor_full = rotor; %used to call get_lambda_up
rotor = rotor(2); %in case the generic rotor struct is given
sigma = rotor.solidity;
cl_a = rotor.aero.cl_alpha;
inner = r<rotor.rd; %inside the downwash circle condition (to select elements in array later)
outer = r>=rotor.rd; %outside the downwash circle condition (to select elements in array later)
%idx_inner = find(r<rotor.rd);

lambda_tot_up_inner = get_lambda_up(F(inner),r(inner),pitch(inner),rotor_full,flowfield);
lambda_up_inner = lambda_tot_up_inner-lambda_inf;
lambda_bot_inner = sqrt((sigma*cl_a*1./(16*F(inner))-lambda_up_inner/(2*rotor.rd^2)-lambda_inf/2).^2+...
    sigma*cl_a*pitch(inner).*r(inner)*1./(8*F(inner)))-sigma*cl_a*1./(16*F(inner))+lambda_up_inner/(2*rotor.rd^2)+lambda_inf/2;

lambda_bot_outer = sqrt((sigma*cl_a*1./(16*F(outer))-lambda_inf/2).^2+sigma*cl_a*pitch(outer).*r(outer)*1./(8*F(outer)))-sigma*cl_a*1./(16*F(outer))+lambda_inf/2;

lambda_bot = [lambda_bot_inner lambda_bot_outer];

end