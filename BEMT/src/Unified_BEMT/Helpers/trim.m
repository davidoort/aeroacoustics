function [collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_or_pitch,trimvar,method,timelimit)
%{
TRIM THE COAXIAL SYSTEM TO A DESIRED CT or UPPER COLLECTIVE SETTING

There are two trim routines in this function: 

- "yaw" simply trims the coaxial around a certain collective angle of both
rotors. Probably the collective of both rotors will be around the average
of the collective_u and collective_l passed initially.

- "thrust" this trims the coaxial rotors to a specific thrust output 
(at torque balance)

This function is used extensively to create validation plots and I have run
into many issues caused by the instability of this fixed point method with
a tweaked proportionality parameter. This caused issues when trimming in
different flight conditions (like high axial velocities) or when using
different methods (leishman vs airfoil).

The rewrite of the code in August was done with the following major
differences:

- Instead of a constant proportionality parameter k I use a
condition-,method-specific gradient calculated in getGradient.

- Add CT, CP (in case you ever need to trim to a specific CP) and
net_torque_coeff to coaxial.state so that getGradient can access those
quantities (+collective_u and collective_l) and calculate what happens if
you input collective_u+dcol, collective_l+dcol or (if the gradient we are 
looking for is with respect to a pitch differential)
collective_u+dcol,collective_l-dcol

Inputs:
    coaxial - (struct object) with operational (state) and geometric
    variables of the coaxial rotor

    atm - (struct object) with atmospheric parameters such as air density

    epsilon - (scalar) convergence accuracy for F-lambda iteration

    CT_or_pitch - (scalar) can be either collective_u (upper rotor pitch 
    angle in degrees) or the desired thrust coefficient of the coaxial
    system

    trimvar - (string) indicates if CT desired or collective desired is being 
    specified

Outputs:
    collective_u - (scalar) pitch angle in degrees of the upper rotor

    collective_l - (scalar) pitch angle in degrees of the lower rotor

    net_torque_dimensional - (scalar) net torque in Nm of the coaxial
    system

    CT - (scalar) non-dimensional thrust coefficient of the coaxial rotor

Other m-files required: 

    BEMT

    trim_torque

    getGradient

MAT-files required: none

Literature referenced: none

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 15-August-2019
%}

%------------- BEGIN CODE --------------

tic %begin trim time-counting, no matter the routine (Torque or CT&Torque)

debug = false;
plots = false;
verbose = false;

%Delta keep in mind that this is currently a lower bound on the exploration
%collective angle. I think delta should be a function of the current
%collective angle.
delta = 0.1; %[deg] what collective angle change is used to calculate the gradient - the smaller, the more local the derivative 

if strcmpi(trimvar,"yaw")
    %do a trim procedure for airfoil method as well - use the delta method
    %of the previous iteration (delta_pitch/delta_thrust instead of k).
    %This theoretically also solves the problem of decreasing pitch to
    %decrease torque (which for negative angles of attack will just spiral
    %down).
    
    coaxial.state.collective_u = CT_or_pitch; %deg
    coaxial.state.collective_l = CT_or_pitch; %deg
    
    %coaxial.state.collective = collective_u; %deg
    [net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method,timelimit,delta); %deg,Nm,-%           (coaxial,atm,epsilon,method,timelimit)



elseif strcmpi(trimvar,"thrust") 
    
    if CT_or_pitch > 1
        error([num2str(CT_or_pitch), ' is an unrealistic value for CT'])
    end

    k1 = 10;
    eps1 = 0.00001;
        
    CT_des = CT_or_pitch;

    %% Begin iteration
    
    [net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method,timelimit,delta); %deg,Nm,-
    

    thrust_error = CT_des-sum(coaxial.state.CT);
     
    
    while toc<timelimit && abs(thrust_error) > eps1
        if strcmpi(method,'vitleish')
            correction = k1*thrust_error;
            delta = 0;
        else
            grad_thrust = getGradient('thrust',coaxial,atm,epsilon,method,debug,delta);
            
            correction = thrust_error/grad_thrust;
        end
        
        coaxial.state.collective_u = abs(coaxial.state.collective_u + correction);
        coaxial.state.collective_l = abs(coaxial.state.collective_l + correction);
        
        
        
         if strcmpi(method, 'leishman') || strcmpi(method, 'vitleish')
            %this might not completely solve the problem with twisted blades
            %and Leishman
            %if coaxial.state.collective_u < 0 || coaxial.state.collective_u < 0
            if coaxial.state.collective_u < delta
                collective_adjustment = delta-coaxial.state.collective_l;
                coaxial.state.collective_u = delta;
                
                coaxial.state.collective_l = coaxial.state.collective_l+collective_adjustment;
                
                
            end
            if coaxial.state.collective_l < delta
                collective_adjustment = delta-coaxial.state.collective_l;
                coaxial.state.collective_l = delta;
                
                coaxial.state.collective_u = coaxial.state.collective_u+collective_adjustment;
            end
            
            %Sanity check that BEMT can be run, otherwise a message
            %pops up and collective is increased-not very elegant atm (it doesn't really do anything)
            error_var = true;
            while error_var
                try
                    %tic
                    [~, ~, Moments, CT_BEMT, CP_BEMT, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
                    %toc
                    error_var = false;
                catch e
                    disp(e.message)
                    coaxial.state.collective_u = coaxial.state.collective_u+0.1;
                    coaxial.state.collective_l = coaxial.state.collective_u+0.1;
                end
            end
            
        end
      
        [net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method,timelimit,delta); %deg,Nm,-

        thrust_error = CT_des - sum(coaxial.state.CT);
    end
    
else
    
    error('Please specify either "yaw" or "thrust" as trim conditions')


end

collective_u = coaxial.state.collective_u;
collective_l = coaxial.state.collective_l;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method,timelimit,delta)

%{
TRIM THE LOWER ROTOR TO ACHIEVE TORQUE BALANCE

This function will basically update the trim setting of the coaxial system,
and since it is given a collective_u, it will return the required
collective_l to balance the torque

Inputs:
    coaxial - (struct object) with operational (state) and geometric
    variables of the coaxial rotor

    atm - (struct object) with atmospheric parameters such as air density

    epsilon - (scalar) convergence accuracy for F-lambda iteration

Outputs:
    collective_l - (scalar) pitch angle in degrees of the lower rotor

    net_torque_dimensional - (scalar) net torque in Nm of the coaxial
    system

    CT - (scalar) non-dimensional thrust coefficient of the coaxial rotor

Other m-files required: 

    BEMT_axial

MAT-files required: none

Literature referenced: none

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 2-June-2019
%}

%------------- BEGIN CODE --------------    
    
    %% Tweak for fast & stable convergence    

    k = 10; %proportionality constant k that seems ideal
    eps = 0.1;
    plots = false;
    debug = false;
    verbose = false;
    
    %% Begin iteration    
    
    [~, ~, Moments, CT_BEMT, CP_BEMT, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
    
    coaxial.state.CT = CT_BEMT;
    coaxial.state.CP = CP_BEMT;  
    coaxial.state.net_torque_coeff = net_torque_coeff;
    
    net_torque_dimensional = net_torque_coeff*sum(abs(Moments(:,3)));
    
    
    while toc<timelimit && abs(net_torque_dimensional)>eps
        
        if strcmpi(method,'vitleish')
            correction = k*coaxial.state.net_torque_coeff;
            delta = 0;
        else
            grad_torque = getGradient('net_torque',coaxial,atm,epsilon,method,debug,delta);
            
            correction = coaxial.state.net_torque_coeff/grad_torque;
        end
        
        coaxial.state.collective_u = coaxial.state.collective_u - correction;
        coaxial.state.collective_l = coaxial.state.collective_l + correction;
        
        if strcmpi(method, 'leishman') || strcmpi(method, 'vitleish')
            %this might not completely solve the problem with twisted blades
            %and Leishman
            %if coaxial.state.collective_u < 0 || coaxial.state.collective_u < 0
            if coaxial.state.collective_u < delta
                collective_adjustment = delta-coaxial.state.collective_l;
                coaxial.state.collective_u = delta;
                
                coaxial.state.collective_l = coaxial.state.collective_l+collective_adjustment;
                
                
            end
            if coaxial.state.collective_l < delta
                collective_adjustment = delta-coaxial.state.collective_l;
                coaxial.state.collective_l = delta;
                
                coaxial.state.collective_u = coaxial.state.collective_u+collective_adjustment;
            end
            
            error_var = true;
            while error_var
                try
                    %tic
                    [~, ~, Moments, CT_BEMT, CP_BEMT, net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
                    %toc
                    error_var = false;
                catch e
                    disp(e.message)
                    coaxial.state.collective_u = coaxial.state.collective_u+0.1;
                    coaxial.state.collective_l = coaxial.state.collective_u+0.1;
                end
            end
            
        end
        
        coaxial.state.CT = CT_BEMT;
        coaxial.state.CP = CP_BEMT; 
        coaxial.state.net_torque_coeff = net_torque_coeff;

        net_torque_dimensional = net_torque_coeff*sum(abs(Moments(:,3)));
        
        
%         if old_net_torque_dimensional + net_torque_dimensional < eps
%             warning("Bouncing around between +- net torques. Nudging k...")
%             if k>0.1
%                 k = 0.9*k;
%             else
%                 k = 1.1*k; %it was getting stuck in infinite loops
%             end
%             
%         end
%         if old_net_torque_dimensional - net_torque_dimensional < eps
%             eps = 1.01*eps;
%             warning(['Increasing eps to ', num2str(eps), ' Nm']);
%         end

    end

    CT = coaxial.state.CT;
    
end   


