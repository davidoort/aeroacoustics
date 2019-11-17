classdef Atmosphere < dynamicprops
    %PARS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T
        gamma_air
        R_air
        rho
        c0
        rel_humidity
        p
        kin_visc
    end
    
    methods
        function obj = Atmosphere()
            
            obj.rho = 1.225; %kg/m^3 sea-level air density
            obj.gamma_air = 1.4; %[-] for air (diatomic gas value)
            obj.R_air = 8.314510/0.0289645; %J/kg/K= N m / N s^2 m^-1 /K = m^2/s^2/K
            obj.T = 20; %degrees Celsius 
            obj.c0 = sqrt(obj.gamma_air*obj.R_air*(obj.T+273.15)); %m/s - speed of sound
            obj.rel_humidity = 0.8; %Humidity of 80% is typical for LA
            obj.p = 1.01325e5; %Ambient pressure in Pa
            kin_visc_arr = [1.2462e-5 1.3324e-5 1.4207e-5 1.5111e-5]; %m^2/s - these are for sea level air density I assume
            temp_arr = [-10 0 10 20]; %deg Celsius
            obj.kin_visc = interp1(temp_arr,kin_visc_arr,obj.T); %dynamic viscosity in m^2/s. This comes from airfoiltools. They give different values at different temperatures
            
            %% FIXED
%             obj.T    = 271+15; % K | Temperature
%             obj.rho  = 1.225; % km/m^3 | density
%             obj.a    = sqrt(1.4*8.314510*(obj.T)/0.0289645); % m/s | speed of sound
%             obj.wind = 0;   % m/s | local wind speed
        end
    end
end
