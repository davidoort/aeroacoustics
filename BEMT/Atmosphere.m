classdef Atmosphere < dynamicprops
    %PARS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Temp
        gamma_air
        R_air
        rho
        c0
    end
    
    methods
        function obj = Atmosphere()
            
            obj.rho = 1.225; %kg/m^3 sea-level air density
            obj.gamma_air = 1.4; %[-] for air (diatomic gas value)
            obj.R_air = 8.314510/0.0289645; %J/kg/K= N m / N s^2 m^-1 /K = m^2/s^2/K
            obj.Temp = 20; %degrees Celsius 
            obj.c0 = sqrt(obj.gamma_air*obj.R_air*(obj.Temp+273.15)); %m/s - speed of sound

            
            %% FIXED
%             obj.T    = 271+15; % K | Temperature
%             obj.rho  = 1.225; % km/m^3 | density
%             obj.a    = sqrt(1.4*8.314510*(obj.T)/0.0289645); % m/s | speed of sound
%             obj.wind = 0;   % m/s | local wind speed
        end
    end
end
