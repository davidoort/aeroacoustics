classdef Rotor < dynamicprops
    %PARS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        params
        rotor
        state
        type
    end
    
    methods
        function obj = Rotor()
            obj.name = "Hermes"; 
            obj.type = "coaxial"; % "single" or "coaxial"
            
            obj.state.trim = 1; %1 means that both rotors have the same geometrical pitch, so same collective setting >1 increases pitch of lower wrt to upper
            obj.state.collective = 8; %deg
            obj.state.axial_vel = 0;
            obj.state.forward_vel = 0;
            obj.state.side_vel = 0;
            obj.state.airspeed = @() norm([obj.state.axial_vel,obj.state.forward_vel,obj.state.side_vel]); %m/s 
            obj.state.incidence_deg = @() rad2deg(atan(obj.state.forward_vel/obj.state.axial_vel)); %deg - positive downward
            obj.state.sideslip = @() atan2(obj.state.side_vel,obj.state.forward_vel); %rad - atan2 because otherwise there is a NaN singularity when in hover
            
            %Cyclic state (input) - for now assumes that the cyclic is
            %coupled between the two rotors
            %these are real coefficients that get multiplied by sine and
            %cosine terms which get added to the pitch matrix (a function
            %of r due to blade design and collective and a function of psi due to cyclic input)
            obj.state.cyclic_s = 0;
            obj.state.cyclic_c = 0;
            %% Create the rotors
            if obj.name == "Hermes"
                
                %Talaria's Hermes II rotor
                
                
                
                %% General params
                
                
                obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                obj.params.interrotor_spacing = 0.15;%percentage of rotor radius
                
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades
                obj.rotor(1).R = 1.28;             %m - radius of blade
                obj.rotor(1).omega = 1200*2*pi/60;
                obj.rotor(1).chord = 0.14;             %m - chord
                obj.rotor(1).solidity = obj.rotor(1).Nb*obj.rotor(1).chord/(pi*obj.rotor(1).R);   % blade solidity - untapered
                obj.rotor(1).hub_radial_fraction = 0.15;
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
                obj.state.collective = 10; % [deg] collective
                obj.rotor(1).pitch_root = 20; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'linear'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.212;         %1/rad - Lift slope, NACA 23015
                obj.rotor(1).aero.alpha_0 = 0; %rad - guess!
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                %% Bottom rotor
                obj.rotor(2) = obj.rotor(1);
                
                
                
                
                
            end
            
            
           
           
        end
        
        
      
        
    end
end
