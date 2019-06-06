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
            obj.name = "Harrington1"; % can be Harrington1, Harrington2, NACA_single, NACA_coax, Hermes
            obj.type = "coaxial"; % "single" or "coaxial"
            
            obj.state.trim = 1; %1 means that both rotors have the same geometrical pitch, so same collective setting >1 increases pitch of lower wrt to upper
            obj.state.pitchdeg = 8; %deg
            obj.state.axial_vel = 0;
            obj.state.tangent_vel = 0;
            obj.state.airspeed = @() norm([obj.state.axial_vel,obj.state.tangent_vel]); %m/s 
            obj.state.incidence_deg = @() rad2deg(atan(obj.state.tangent_vel/obj.state.axial_vel)); %deg - positive downward
            
            %% Create the rotors
            if obj.name == "Harrington1"
                
                
                %% General params
            

                %obj.type = "single"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 1 had an interrotor spacing of 0.186 R = 1.41 m
            
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(1).R = 7.62;             %m - arbitrary
                obj.rotor(1).omega = 15;            %rad/s
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.027;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(2).R = 7.62;             %m 
                obj.rotor(2).omega = 15; %rad/s - arbitrary
                obj.rotor(2).chord = 0.15;             %m - arbitrary
                obj.rotor(2).solidity = 0.027;
                
                
                % Aerodynamics
                obj.rotor(2).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
              
                
            elseif obj.name == "Harrington2"
                
                
                %% General params
            

                %obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m 
            
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(1).R = 7.62;             %m
                obj.rotor(1).omega = 15;        %rad/s
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.076;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(2).R = 7.62;             %m - arbitrary
                obj.rotor(2).omega = 15;        %rad/s - arbitrary
                obj.rotor(2).chord = 0.15;             %m - arbitrary
                obj.rotor(2).solidity = 0.076;
                
                
                % Aerodynamics
                obj.rotor(2).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
                
            elseif obj.name == "NACA_single"
                
                
                %% General params
            

                obj.type = "single"; % "single" or "coaxial"
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m 
            
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 3;             % Number of blades - arbitrary
                obj.rotor(1).R = 1.4859;             %m
                obj.rotor(1).omega = 249.022/obj.rotor(1).R*2*pi; %the magic number is the tip speed in m/s from Leishman
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.2292;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5;         %1/rad - Lift slope, NACA 16-series... GUESS at the moment
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
       
                
            elseif obj.name == "NACA_coax"
                
                
                %% General params
            

                obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                %Rotor 2 had an interrotor spacing of 0.16 R = 1.21 m 
            
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 3;             % Number of blades - arbitrary
                obj.rotor(1).R = 1.524;             %m
                obj.rotor(1).omega = 249.022/obj.rotor(1).R*2*pi; %the magic number is the tip speed in m/s from Leishman
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.078;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.11;         %1/rad - Lift slope, Clark-Y airfoil GUESS!
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2).Nb = 3;             % Number of blades - arbitrary
                obj.rotor(2).R = 7.62;             %m - arbitrary
                obj.rotor(1).omega = 249.022/obj.rotor(1).R*2*pi; %the magic number is the tip speed in m/s from Leishman
                obj.rotor(2).chord = 0.15;             %m - arbitrary
                obj.rotor(2).solidity = 0.078;
                
                
                % Aerodynamics
                obj.rotor(2).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
                
            elseif obj.name == "Hermes"
                
                %Talaria's Hermes II rotor
                
                
                  
                %% General params
            

                obj.type = "coaxial"; % "single" or "coaxial"
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
            
            
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades
                obj.rotor(1).R = 1.28;             %m - radius of blade
                obj.rotor(1).omega = 1200*2*pi/60;
                obj.rotor(1).chord = 0.14;             %m - chord
                obj.rotor(1).solidity = obj.rotor(1).Nb*obj.rotor(1).chord/(pi*obj.rotor(1).R);   % blade solidity - untapered
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.212;         %1/rad - Lift slope, NACA 23015
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
               
                
                %% Bottom rotor
                obj.rotor(2).Nb = 2;  %# Number of blades
                obj.rotor(2).R = 1.28;  %m - using reference values from papers for now
                obj.rotor(2).omega = 1200*2*pi/60;
                obj.rotor(2).chord = 0.14;  %m - using reference values from papers for now
                obj.rotor(2).solidity = obj.rotor(2).Nb*obj.rotor(2).chord/(pi*obj.rotor(2).R);   % aero solidity - not sure what to take as chord for tapered aero
                
                
                % Aerodynamics - maybe a bit of overkill at the moment
                obj.rotor(2).aero.cl_alpha = 5.212;         %1/rad - Lift slope, NACA 23015
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
             
                
                
            end
            
            
           
           
        end
        
        
      
        
    end
end
