classdef Observer < dynamicprops
    %Properties of the observer such as relative position and relative
    %motion wrt to the rotor, important for acoustic computations.
    
    %DISTANCES AND ANGLES WRT TO HUB OF LOWER ROTOR
    
    properties
        r1
        r
        x
        theta
        psi
    end
    
    methods
        function obj = Observer()
            
            %I need a nicer way of defining the observer position. The
            %rotor is at coordinate (0,0,0) pointing in positive x-direction in the rotor body frame (along the rotational axis)
            
            
            obj.r1      = 10; % [m] - see the definition in Hanson's paper "Noise of Counter-Rotation Propellers"
            obj.theta   = 90; % [deg] - %defined as the angle between the rotation axis of the rotor and the line connecting the lower rotor hub with the observer as pictorially depicted in Fig.3
            obj.x       = 10; % [m] - see the definition in Hanson's paper "Noise of Counter-Rotation Propellers"
            obj.r       = @() norm([obj.x,obj.r1]); % [m] - see the definition in Hanson's paper "Noise of Counter-Rotation Propellers"
            obj.psi     =  0;   %[rad]
            %from theta and axial and tangent velocity it will be possible
            %to calculate M_x and M_T
        end
    end
end
