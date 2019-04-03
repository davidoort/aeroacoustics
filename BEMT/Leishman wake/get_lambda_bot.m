function lambda_bot = get_lambda_bot(F,r,pitch,rotor,flowfield,lambda_u)

%F = Prandtl tip loss (array)
%r = radial position of blade (array)
%pitch of blade element wrt to rotor plane in rad (array)
%rotor = rotor struct

lambda_inf = flowfield(2).lambda_inf;

rotor = rotor(2); %in case the generic rotor struct is given
sigma = rotor.solidity;
cl_a = rotor.aero.cl_alpha;
inner = r<rotor.rd; %inside the downwash circle condition (to select elements in array later)
outer = r>=rotor.rd; %outside the downwash circle condition (to select elements in array later)
%idx_inner = find(r<rotor.rd);

%{
Make sure that lambda bottom inner uses lambda_u from root to tip instead
of from root to rd - "Note that the contracting streamtubes must be mapped
from originating points on the upper rotor to receiving points on the
lower rotor within the inner area that is affected by the wake, i.e., a
streamtube emanating from radial point r on the upper rotor will map to
new radial point ar on the lower rotor"
%}

%Procedure has been verified in the command window. lambda_u now has the
%same shape but with a coarser resolution to match the dimensions of
%r(inner).
interp_length = length(r(inner)); 
ri = linspace(0,1,interp_length);
lambda_u = interp1(r,lambda_u,ri);

%since lambda inf is uniform, it can just be truncated so lambda_inf(inner)

lambda_bot_inner = sqrt((sigma*cl_a*1./(16*F(inner))-lambda_u/(2*rotor.rd^2)-lambda_inf(inner)/2).^2+...
    sigma*cl_a*pitch(inner).*r(inner)*1./(8*F(inner)))-sigma*cl_a*1./(16*F(inner))+lambda_u/(2*rotor.rd^2)+lambda_inf(inner)/2;

lambda_bot_outer = sqrt((sigma*cl_a*1./(16*F(outer))-lambda_inf(outer)/2).^2+sigma*cl_a*pitch(outer).*r(outer)*1./(8*F(outer)))-sigma*cl_a*1./(16*F(outer))+lambda_inf(outer)/2;

lambda_bot = [lambda_bot_inner lambda_bot_outer];

end