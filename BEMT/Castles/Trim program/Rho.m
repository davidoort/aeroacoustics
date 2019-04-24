% Calculates the density of air at a specific altitude

function density = rho(h)

% Define globals
global hstrat rho0 lambda T0 g R0 h

if h < hstrat
	density   = rho0*(1+lambda*h/T0)^(-1*(g/(R0*lambda)+1));
else
	Ts        = T0+lambda*hstrat;
	rho_strat = rho0*(1+lambda*hstrat/T0)^(-1*(g/(R0*lambda)+1));
	density   = rho_strat*exp(-g/R0/Ts*(h-hstrat));
end