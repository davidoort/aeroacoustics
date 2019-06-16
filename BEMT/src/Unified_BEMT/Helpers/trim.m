function [collective_u, collective_l, net_torque_dimensional, CT] = trim(coaxial,atm,epsilon,CT_or_pitch,trimvar,method)
%{
TRIM THE COAXIAL SYSTEM TO A DESIRED CT or UPPER COLLECTIVE SETTING

The inputs of this function change depending on whether you want to trim
the rotor at a specified value of CT or at a specified upper rotor pitch
angle. This could easily be extended to take lower rotor pitch angle as
well but it is not necessary for the rest of the pipeline at the moment.

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

    BEMT_axial

    trim_torque

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
debug = false;

if strcmpi(trimvar,"pitch_upper")
    
    collective_u = CT_or_pitch; %deg
    coaxial.state.collective = collective_u; %deg
    [collective_l,net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method); %deg,Nm,-
    
elseif strcmpi(trimvar,"CT")
    
    if CT_or_pitch > 1
        error([num2str(CT_or_pitch), ' is an unrealistic value for CT'])
    end

    %% Tweak for fast & stable convergence 

    k1 = 70;
    eps1 = 0.00001;
    plots = false;
    verbose = false;
        
    CT_des = CT_or_pitch;

    %% Begin iteration
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
    
    thrust_error = CT_des - coaxial.state.CT;
    
    while abs(thrust_error) > eps1
        
        coaxial.state.collective = coaxial.state.collective + k1*thrust_error;
        
        if coaxial.state.collective<0
            disp('Collective set to 0 instead of negative')
            coaxial.state.collective = 0;
            
        end
        
        
        [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
            coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
        
        [collective_l,net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method); %deg,Nm,-
        
        thrust_error = CT_des - coaxial.state.CT;
        
    end
    
    collective_u = coaxial.state.collective;
    
else
    error('Please specify either "CT" or "pitch_upper" as trim conditions')


end


function [collective_l,net_torque_dimensional,CT] = trim_torque(coaxial,atm,epsilon,method)

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

    k = 1; %proportionality constant k that seems ideal - what also worked was doing k*coaxial.state.net_torque_coeff even if I don't change epsilon
    eps = 0.1;
    plots = false;
    verbose = false;
    
    %% Begin iteration    
    
    [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
    
    net_torque_dimensional = coaxial.state.net_torque_coeff*coaxial.state.torque;
    
    while abs(net_torque_dimensional)>eps
        old_net_torque_dimensional = net_torque_dimensional;
        
        coaxial.state.trim = coaxial.state.trim + k*coaxial.state.net_torque_coeff;
        if coaxial.state.trim < 0
           warning(strcat("Trim value = ", num2str(coaxial.state.trim), " resetting trim to 1 and lowering k"))
           k = k*0.9;
        end
        
        %tic
        [coaxial.state.thrust, coaxial.state.torque, coaxial.state.power, ...
        coaxial.state.CT, coaxial.state.CP, coaxial.state.net_torque_coeff] = BEMT(coaxial,atm,epsilon,plots,verbose,method,debug);
        %toc
        net_torque_dimensional = coaxial.state.net_torque_coeff*coaxial.state.torque;
        
        if abs(old_net_torque_dimensional) - abs(net_torque_dimensional) < eps
            warning("Bouncing around between +- net torques. Nudging k...")
            k = 1.01*k; %it was getting stuck in infinite loops 
        end
    end
    
    collective_l = coaxial.state.trim*coaxial.state.collective;
    CT = coaxial.state.CT;
    
end   

end
