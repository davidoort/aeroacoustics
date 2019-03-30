function Fcf= Prandtl_tip_loss(r,lambda,rotor)
%Check of they neatly wrote inputs and outputs in robot dynamics course

%r is the radial position of the blade element, lambda the inflow ratio at
%the element. The third argument should either be rotor(1) or rotor(2)

Fcf = (2/pi)*acos(exp(-rotor.Nb*(1-r)/(2*lambda)));

end
