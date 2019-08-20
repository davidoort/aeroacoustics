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
            obj.name = "Hermes"; % can be Harrington1, Harrington2, NACA_single, NACA_coax, Hermes, Bumblebee_wing, Bumblebee_canard
            obj.type = "coaxial"; % "single" or "coaxial"
            
            obj.state.trim = 1; %1 means that both rotors have the same geometrical pitch, so same collective setting >1 increases pitch of lower wrt to upper
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
            if obj.name == "Harrington1"
                
                
                %% General params
            

                %obj.type = "single"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 1 had an interrotor spacing of 0.186 R = 1.41 m
                %IF LANDGREBE (actually I need it anyway):
                obj.params.interrotor_spacing = 0.186; % [-] fraction of upper rotor radius
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - from Harrington paper
                obj.rotor(1).R = 3.81;             %m - from Harrington paper (12.5 ft)
                obj.rotor(1).omega = 500/12.5;            %rad/s - Fig 7 harrington paper
                obj.rotor(1).root_chord = 0.28702;             %m  - from Harrington paper
                obj.rotor(1).tip_chord = 0.11176;             %m - from Harrington paper
                obj.rotor(1).solidity = 0.027; %from Harrington paper
                obj.rotor(1).hub_radial_fraction = 20/150; %from Harrington paper
                
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 twist
               
                obj.rotor(1).pitch_root = 10; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'linear'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                %obj.rotor(1).airfoil.name = "NACA0012"; %airfoil instead
                %of aero...
                obj.rotor(1).airfoil.name = "NACA0012";
                obj.rotor(1).aero.cl_alpha = 5.73;         %1/rad - Lift slope, NACA 0012 from Harrington report
                %obj.rotor(1).aero.cl_alpha = 2*pi; doesn't make much of a
                %difference to verification/validation plots
                obj.rotor(1).aero.alpha_0 = 0; %rad
                obj.rotor(1).aero.Cd0 =  0.007; %0.011
                obj.rotor(1).aero.D1 = 0; %1/rad
                obj.rotor(1).aero.D2 = 1; %1/rad^2
                
                
                
                %% Bottom rotor
                
                obj.rotor(2) = obj.rotor(1);
              
                
            elseif obj.name == "Harrington2"
                
                
                %% General params
            

                %obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m 
                obj.params.interrotor_spacing = 0.16; % [-] fraction of upper rotor radius
                
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - from Harrington paper
                obj.rotor(1).R = 3.81;             %m - from Harrington paper (12.5 ft)
                obj.rotor(1).omega = (327+392)/2/12.5;        %rad/s - Fig 8 harrington paper
                obj.rotor(1).root_chord = 0.4572;             %m  - from Harrington paper
                obj.rotor(1).tip_chord = 0.4572;             %m  - from Harrington paper
                obj.rotor(1).solidity = 0.076; %from Harrington paper
                obj.rotor(1).hub_radial_fraction = 30/150; %from Harrington paper
                 
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
              
                obj.rotor(1).pitch_root = 20; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'linear'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                obj.rotor(1).airfoil.name = "NACA0012";
                obj.rotor(1).aero.cl_alpha = 5.73;         %1/rad - Lift slope, NACA 0012 from Harrington report
                obj.rotor(1).aero.alpha_0 = 0; %rad
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 1;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2) = obj.rotor(1);
                
            elseif obj.name == "NACA_single"
                %Prop information
                %{
                
                "This NACA propeller has thin blades (about 8% thickness-to-chord ratio)".
                "The three-blade propeller used for this investigation had ~ diameter of
                9.75 feet. The blades were made from solid 6415 steel nnd had SYMMETRICAL
                16-series airfoil sections (will use this airfoil data
                http://airfoiltools.com/airfoil/details?airfoil=naca16006-il). - the
                Reynolds number I have to figure out by looking at the effective reynolds
                number at each advance ratio...
                The blade-thickness ratio varied from 0.06 at the spinner to 0,02 at the tip,
                and the chord was constant for the length of the exposed blade."

                No information about pitch distribution except in a different paper, where
                it mentiones that "The pitch distribution was the same for the three
                propellers and was such that the chord line of each section lay in a plane
                parallel to the resultant section velocity corrected for the induced flow
                for the design advance ratio, J = 2.0. The induced flow was calculated
                according to Goldstein's theory for straight propellers." Am I supposed to
                implement Goldstein's theory to find the pitch distribution?! Closest thing
                I could do now is have a pitch distribution that is the same as alpha when
                I spin it with untwisted blades with a forward flight ratio of 2.0/pi
                (J/pi). Actually, that gives a non-sensical twist...

                                v-tip = 249.022;
                
                %}

                
                %% General params
            

                obj.type = "single"; % by definition
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                
           
                %% Top rotor
                
                obj.rotor(1).Nb = 3;             % Number of blades 
                obj.rotor(1).R = 1.4859;             %m
                obj.rotor(1).omega = 249.022/obj.rotor(1).R; %the magic number is the tip speed in m/s from Leishman
                obj.rotor(1).solidity = 0.2292;
                obj.rotor(1).chord = pi*obj.rotor(1).solidity*obj.rotor(1).R/obj.rotor(1).Nb; %m - calculated, since we know it is constant and = sigma*pi*R/B
                obj.rotor(1).hub_radial_fraction = 0.1; %guess
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 twist
                
                obj.rotor(1).pitch_root = 50; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'linear'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 6.7;         %1/rad - Lift slope, NACA 16-series... GUESS at the moment
                obj.rotor(1).aero.alpha_0 = 0; %rad - symmetrical
                obj.rotor(1).aero.Cd0 =  0.003; %0.003
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 3; %this could be causing problems? 
                %very high profile drag apparently in the axial validation
                %plots. When I run the code at axial vel of 210 and I try
                %to get a CT = 0.00037 then it trims to a collective of 74.
                %In the disks I can see there is negative thrust around the
                %center of the rotor... Found it!! Wrong choice of equation
                %for CT!!
       
                
            elseif obj.name == "NACA_coax"
                
                
                %% General params
                obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                obj.params.interrotor_spacing = 0.8585/1.524; %fraction of the radius from Bierman
            

                %% Top rotor
                obj.rotor(1).Nb = 3;             % Number of blades - from the picture in Bierman
                obj.rotor(1).R = 1.524;             %m
                obj.rotor(1).omega = 87.47/obj.rotor(1).R; %the magic number is the tip speed in m/s from Leishman
                obj.rotor(1).solidity = 0.078; %try with this solidity and double this one
                obj.rotor(1).chord = pi*obj.rotor(1).solidity*obj.rotor(1).R/obj.rotor(1).Nb; %m - average chord (not constant this time); chord = sigma*pi*R/B
                obj.rotor(1).hub_radial_fraction = 0.2; %from figure 1 in Bierman
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
  
                obj.rotor(1).pitch_root = 20; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'linear'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics - HS (Hamilton Standard) 3155-6 and HS 3156-6 blades
                obj.rotor(1).aero.cl_alpha = 6.1;         %1/rad - Lift slope, Clark-Y airfoil!
                obj.rotor(1).aero.alpha_0 = deg2rad(-3.75); %rad - from airfoil data
                obj.rotor(1).aero.Cd0 =  0.007;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.8;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2) = obj.rotor(1);
                
            elseif obj.name == "Hermes"
                
                %Talaria's Hermes II rotor
                
                
                  
                %% General params
            

                obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                obj.params.interrotor_spacing = 0.4; 

            
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades
                obj.rotor(1).R = 1.28;             %m - radius of blade
                obj.rotor(1).omega = 1200*2*pi/60;
                obj.rotor(1).chord = 0.14;             %m - chord
                obj.rotor(1).root_chord = 0.14;             %m - chord
                obj.rotor(1).tip_chord = 0.14;             %m - chord
                obj.rotor(1).solidity = obj.rotor(1).Nb*obj.rotor(1).chord/(pi*obj.rotor(1).R);   % blade solidity - untapered
                obj.rotor(1).hub_radial_fraction = 0.15;
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
             
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
                
            elseif obj.name == "Bumblebee_wing"
                
                
                %% General params
            
                
                obj.type = "single"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m 
            
                
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - from Harrington paper
                obj.rotor(1).R = 1;             %m - from Harrington paper (12.5 ft)
                obj.rotor(1).omega = 300;        %rad/s - Fig 8 harrington paper
                obj.rotor(1).root_chord = 0.1;             %m  - from Harrington paper
                obj.rotor(1).tip_chord = 0.1;             %m  - from Harrington paper
                obj.rotor(1).solidity = 0.16; %from Harrington paper
                obj.rotor(1).hub_radial_fraction = 0.15; %from Harrington paper
                 
                
                 % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
     
                obj.rotor(1).pitch_root = 14; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'ideal'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                %obj.rotor(1).airfoil.name = "NACA0012"; %airfoil instead
                %of aero...
                obj.rotor(1).airfoil.name = "NACA0012";
                obj.rotor(1).aero.cl_alpha = 5.73;         %1/rad - Lift slope, NACA 0012 from Harrington report
                %obj.rotor(1).aero.cl_alpha = 2*pi; doesn't make much of a
                %difference to verification/validation plots
                obj.rotor(1).aero.alpha_0 = 0; %rad
                obj.rotor(1).aero.Cd0 =  0.007;
                obj.rotor(1).aero.D1 = 0;
                obj.rotor(1).aero.D2 = 1;
                
            elseif obj.name == "Bumblebee_canard"
                
                
                %% General params
                
                
                obj.type = "single"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m
                
                
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - from Harrington paper
                obj.rotor(1).R = 1;             %m - from Harrington paper (12.5 ft)
                obj.rotor(1).omega = 300;        %rad/s - Fig 8 harrington paper
                obj.rotor(1).root_chord = 0.1;             %m  - from Harrington paper
                obj.rotor(1).tip_chord = 0.1;             %m  - from Harrington paper
                obj.rotor(1).solidity = 0.181; %from Harrington paper
                obj.rotor(1).hub_radial_fraction = 0.15; %from Harrington paper
                
                
                % geometric pitch - IF untwisted blade simply write
                % 'linear' with 0 wtist
             
                obj.rotor(1).pitch_root = 16; % [deg] - for an ideal twist it will change the relative angles between root and tip
                obj.rotor(1).twist_type = 'ideal'; % 'linear', 'ideal'
                obj.rotor(1).twistdeg = 0; %[deg] for a fixed collective, this is the same as changing the rate of twist
                
                % Aerodynamics
                %obj.rotor(1).airfoil.name = "NACA0012"; %airfoil instead
                %of aero...
                obj.rotor(1).airfoil.name = "NACA0012";
                obj.rotor(1).aero.cl_alpha = 5.73;         %1/rad - Lift slope, NACA 0012 from Harrington report
                %obj.rotor(1).aero.cl_alpha = 2*pi; doesn't make much of a
                %difference to verification/validation plots
                obj.rotor(1).aero.alpha_0 = 0; %rad
                obj.rotor(1).aero.Cd0 =  0.007;
                obj.rotor(1).aero.D1 = 0;
                obj.rotor(1).aero.D2 = 1;
                
                
                
                
                
            end
            
            
           
           
        end
        
        
      
        
    end
end
