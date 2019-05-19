classdef Coaxial < dynamicprops
    %PARS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        params
        rotor
        state
        
    end
    
    methods
        function obj = Coaxial()
            obj.name = "Hermes"; % can be Harrington, Hermes
            
            
            obj.state.trim = 1; %1 means that both rotors have the same geometrical pitch, so same collective setting >1 increases pitch of lower wrt to upper
            obj.state.pitchdeg = 10; %deg
            obj.state.airspeed = 0; %m/s 
            obj.state.incidence_deg = 0; %deg - positive downward
            obj.state.axial_vel = obj.state.airspeed*cos(deg2rad(obj.state.incidence_deg));
            obj.state.tangent_vel = obj.state.airspeed*sin(deg2rad(obj.state.incidence_deg));
            
            
            %% Create the rotors
            if obj.name == "Harrington"
                
                
                %% General params
            

            
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
        
            
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(1).R = 1.8;             %m - arbitrary
                obj.rotor(1).rpm = 800;             %rev/min - arbitrary
                obj.rotor(1).omega = obj.rotor(1).rpm*2*pi/60;
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.027;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(2).R = 1.8;             %m - arbitrary
                obj.rotor(2).rpm = 800;  %rev/min - arbitrary
                obj.rotor(2).omega = obj.rotor(2).rpm*2*pi/60; %rad/s - arbitrary
                obj.rotor(2).chord = 0.15;             %m - arbitrary
                obj.rotor(2).solidity = 0.027;
                
                
                % Aerodynamics
                obj.rotor(2).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
              
                
                
            elseif obj.name == "Hermes"
                
                %Talaria's Hermes II rotor
                
                
                  
                %% General params
            

            
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
            
            
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades
                obj.rotor(1).R = 1.28;             %m - radius of blade
                obj.rotor(1).rpm = 1200;             %rev/min - angular speed
                obj.rotor(1).omega = obj.rotor(1).rpm*2*pi/60;
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
                obj.rotor(2).rpm = 1200;  %rev/min - using reference values from papers for now
                obj.rotor(2).omega = obj.rotor(2).rpm*2*pi/60;
                obj.rotor(2).chord = 0.14;  %m - using reference values from papers for now
                obj.rotor(2).solidity = obj.rotor(2).Nb*obj.rotor(2).chord/(pi*obj.rotor(2).R);   % aero solidity - not sure what to take as chord for tapered aero
                
                
                % Aerodynamics - maybe a bit of overkill at the moment
                obj.rotor(2).aero.cl_alpha = 5.212;         %1/rad - Lift slope, NACA 23015
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
             
                
                
            end
            
            
           
           
        end
        
        function updateDependenciesPU(obj) % at the moment this will update everything even though some things might not have changed
            
            obj.state.axial_vel = obj.state.airspeed*cos(deg2rad(obj.state.AoAdeg));
            obj.state.tangent_vel = obj.state.airspeed*sin(deg2rad(obj.state.AoAdeg));
            %% Update the rotors
            if obj.name == "Harrington"
                
                
                %% General params
                
                
                
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                
                
                
                %Validation with Harrington Rotor 1 - things assumed are pitch of rotor and
                %CL_alpha
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(1).R = 1.8;             %m - arbitrary
                obj.rotor(1).rpm = 800;             %rev/min - arbitrary
                obj.rotor(1).omega = obj.rotor(1).rpm*2*pi/60;
                obj.rotor(1).chord = 0.15;             %m - arbitrary
                obj.rotor(1).solidity = 0.027;
                
                % Aerodynamics
                obj.rotor(1).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(1).aero.Cd0 =  0.011;
                obj.rotor(1).aero.D1 = 0.0;
                obj.rotor(1).aero.D2 = 0.0;
                
                
                
                %% Bottom rotor
                
                obj.rotor(2).Nb = 2;             % Number of blades - arbitrary
                obj.rotor(2).R = 1.8;             %m - arbitrary
                obj.rotor(2).rpm = 800;  %rev/min - arbitrary
                obj.rotor(2).omega = obj.rotor(2).rpm*2*pi/60; %rad/s - arbitrary
                obj.rotor(2).chord = 0.15;             %m - arbitrary
                obj.rotor(2).solidity = 0.027;
                
                
                % Aerodynamics
                obj.rotor(2).aero.cl_alpha = 5.156;         %1/rad - Lift slope, NACA 0012 from graph
                obj.rotor(2).aero.Cd0 =  0.011;
                obj.rotor(2).aero.D1 = 0.0;
                obj.rotor(2).aero.D2 = 0.0;
                
                
                
            elseif obj.name == "Hermes"
                
                %Talaria's Hermes II rotor
                
                
                
                %% General params
                
                
                
                obj.params.kappaint = 1.28;
                obj.params.kappa = 1.15;
                obj.params.rd = 0.82; %[non-dimensionalised by R] annulus - "assumption consistent with the results obtained by Leishman using the free wake method"
                
                
                
                %% Top rotor
                
                obj.rotor(1).Nb = 2;             % Number of blades
                obj.rotor(1).R = 1.28;             %m - radius of blade
                obj.rotor(1).rpm = 1200;             %rev/min - angular speed
                obj.rotor(1).omega = obj.rotor(1).rpm*2*pi/60;
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
                obj.rotor(2).rpm = 1200;  %rev/min - using reference values from papers for now
                obj.rotor(2).omega = obj.rotor(2).rpm*2*pi/60;
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
