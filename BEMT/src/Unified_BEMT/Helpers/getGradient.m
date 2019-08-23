function grad = getGradient(grad_type,coaxial,atm,epsilon,method,debug)
%{
Preliminary function that gives the numerical derivative of thrust/torque 
wrt collective about the current collective angle

Currently this function evaluates BEMT twice (for simplicity) but it can be
optimized to take the current CT, CP and net_torque and just evaluate BEMT
once

%}

plots = false;
acoustics = false;
verbose = false;
delta = 0.001;


if strcmpi(grad_type,'net_torque')
    
    net_torque_old = coaxial.state.net_torque_coeff; %linearization point CT
    collective_u = coaxial.state.collective_u;
    collective_l = coaxial.state.collective_l;
    
    coaxial.state.collective_u = collective_u+delta; 
    coaxial.state.collective_l = collective_l-delta;
    
    [~, ~, ~, ~, ~, net_torque_new] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);

    
    %gradient is the change in net torque coefficient caused by an increase
    %in the upper collective and a decrease in the lower collective
    
    grad = (net_torque_new-net_torque_old)/delta;
    
    %Reset the state of the coaxial (even though a displacement by delta is not very significant)
    coaxial.state.collective_u = collective_u; 
    coaxial.state.collective_l = collective_l;

    
elseif strcmpi(grad_type,'thrust')
    
    CT_old = coaxial.state.CT; %linearization point CT
    collective_u = coaxial.state.collective_u;
    collective_l = coaxial.state.collective_l;
    
    coaxial.state.collective_u = collective_u+delta; 
    coaxial.state.collective_l = collective_l+delta;
    
    [~, ~, ~, CT_new, ~, ~] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug,acoustics);
    
    
    %gradient is the change in total thrust coefficient (CT_BEMT) caused by an increase in collective
    %of both rotors by a delta amount
    
    grad = (sum(CT_new)-sum(CT_old))/delta;
    
    %Reset the state of the coaxial (even though a displacement by delta is not very significant)
    coaxial.state.collective_u = collective_u; 
    coaxial.state.collective_l = collective_l;

end
