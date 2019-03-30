function Fcf= Prandtl_tip_loss(r,lambda,rotor)
%Check of they neatly wrote inputs and outputs in robot dynamics course

%r is the radial position of the blade element (array)
%lambda the inflow ratio at the element (array) 
%The third argument should either be rotor(1) or rotor(2) 

dim = size(rotor);

if dim(2)>1
    error("Function expects either rotor(1) or rotor(2)")
end

Fcf = (2/pi)*acos(exp(-rotor.Nb*(ones(1,length(r))-r)./(2*lambda)));

end
