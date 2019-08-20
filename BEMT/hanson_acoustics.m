function [dB] = hanson_acoustics(coaxial,atm,r_mat,dz0,CLk,CDk)
%{

Based on "Noise of Counter-Rotation Propellers" by Hanson
This implementation will assume chord compactness and that at t=0 the
blades start aligned with the observer.

dr = dz assumed.

just taking the real part of the expression - assumed. Maybe the norm is
another option (so computing real and imaginary parts)

FOR DEBUGGING PURPOSES, OBSERVER ANGLES ARE IN DEGREES, PSI is still in RAD

Two (probably wrong) approaches have been thought of to deal with non-axial
flight. "angle" violates the definition of angles in the theory of the
paper (by adding an angle) whilst "velocity" violates the definition of velocities Mx and Mt

CLk.upper = CLk_u; ...

%}

p_ref = 2.*10^-5; %Pa

approach = 'angle'; % 'angle' or 'velocity'

% Constants abbreviation - same nomenclature as the paper
rho0 = atm.rho;
c0 = atm.c0;

B(1) = coaxial.rotor(1).Nb;
B(2) = coaxial.rotor(2).Nb;
omega(1) = coaxial.rotor(1).omega;
omega(2)= coaxial.rotor(1).omega;
omega12 = omega(1)/omega(2);
D(1) = 2.*coaxial.rotor(1).R; 
D(2) = 2.*coaxial.rotor(2).R; 
B_D(1) = mean(coaxial.rotor(1).root_chord,coaxial.rotor(1).tip_chord)/D(1);
B_D(2) = mean(coaxial.rotor(2).root_chord,coaxial.rotor(2).tip_chord)/D(2);

observer = Observer();
theta = observer.theta; %defined as the angle between the rotation axis of the rotor and the line connecting the lower rotor hub with the observer as pictorially depicted in Fig.3
r1 = observer.r1;
r = observer.r;

time = linspace(0,1,10); %s


%setup constant matrices

r_vec = r_mat(1,:);
z0 = repmat(r_vec,)

m_range = -10:10; %sound harmonic range
k_range = -10:10; %loading harmonic range - want to cut off frequencies (unless I apply this Gaussian filter)

m = repmat...
k = repmat 


for rotor_idx = [1,2]
    %% Calculate Clk and Cdk
    
    p_t(:,rotor_idx) = zeros(1,length(time));
    


    %Create a time-history for each rotor
    for t_idx = 1:length(time)
        t = time(t_idx);
        if B(1) == B(2) && omega(1) == omega(2)
            %% Special case - acoustic interference? Nonetheless, for verification purposes
            
            %S = coaxial.params.interrotor_spacing.*D(rotor_idx)/2;
            
            %match nomenclature
            omega = omega(1);
            B = B(1);
            
            %% Calculate relative velocities
            
            if strcmpi(approach, 'angle')
                Mx = norm([coaxial.state.axial_vel,coaxial.state.tangent_vel])/c_0;
                Mt = omega.*D(rotor_idx)/2/c_0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta-atan2d(coaxial.state.tangent_vel/coaxial.state.axial_vel); %[deg]
                theta2 = theta; % [deg]
                
            elseif strcmpi(approach, 'velocity')
                Mx = coaxial.state.axial_vel/c_0;
                Mt = (omega.*D(rotor_idx)/2+coaxial.state.tangent_vel)/c_0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta; %[deg]!!
                theta2 = theta; %[deg]!!
            end
            
            %% Equation 6 exponential (for verification)
            
            arg_exp = (m-2.*k).*B(-pi/2)+m.*B(omega.*r/c0-omega.*t); %matrix (k_range x m_range)
            EXP_real = -sin(arg_exp); %matrix
            EXP_imag = cos(arg_exp); %matrix
            
            %% Equation 6 integral (for verification)
            
            % I can do element wise on a three-dimensional matrix (k, m, radius)! Maybe I
            % can add time and rotors to the mix - > 5D!
            kx = 2.*B.*Mt/Mr.*(m/(1-Mx.*cosd(theta1))-2.*k).*B_D(rotor_idx);
            ky = -2.*B./(z0.*Mr).*(m.*(Mr^2.*cosd(theta1)-Mx)/(1-Mx.*cosd(theta1))+2.*k.*Mx).*B_D(rotor_idx); %theta1 is a bit of a guess in the first cosd
            arg_int = Mr^2.*besselj((m-2.*k).*B,m.*B.*z0.*Mt.*sind(theta2)/(1-Mx.*cosd(theta1))).*(kx.*CDk/2+ky.*CLk/2).*dz0; %3D matrix, third dimension being the radius
            INT = sum( arg_int , 3); % 2D matrix
            
        else
            %% Calculate relative velocities
            
            if strcmpi(approach, 'angle')
                Mx = norm([coaxial.state.axial_vel,coaxial.state.tangent_vel])/c_0;
                Mt = omega(rotor_idx).*D(rotor_idx)/2/c_0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta-atan2d(coaxial.state.tangent_vel/coaxial.state.axial_vel); %[deg]
                theta2 = theta; %[deg]
            elseif strcmpi(approach, 'velocity')
                Mx = coaxial.state.axial_vel/c_0;
                Mt = (omega(rotor_idx).*D(rotor_idx)/2+coaxial.state.tangent_vel)/c_0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta; %[deg]
                theta2 = theta; %[deg]
            end
            
            %% Equation 1 exponential (for verification wrt eq.6)
            arg_exp = (m.*B(2)-k.*B(1)).*(-pi/2)+(m.*B(2).*omega(2)+k.*B(1).*omega(1)).*(r/c0-t); %matrix (k_range x m_range)
            EXP_real = -sin(arg_exp); %matrix
            EXP_imag = cos(arg_exp); %matrix
            
            
            %% Equation 1 integral (for verification wrt eq.6)
            
            % I can do element wise on a three-dimensional matrix (k, m, radius)! Maybe I
            % can add time and rotors to the mix - > 5D!
            kx = 2.*Mt/Mr.*((m.*B(2)+k.*B(1).*omega12)/(1-Mx.*cosd(theta1))-k.*B(1).*(1+omega12)).*B_D(rotor_idx);
            ky = -2/Mr.*((m.*B(2)+k.*B(1).*omega12).*Mt^2.*z0.*cosd(theta1)/(1-Mx.*cosd(theta1))-((m.*B(2)-k.*B(1)).*Mx)./z0).*B_D(rotor_idx); %theta1 is a bit of a guess in the first cosd
            arg_int = Mr^2.*besselj((m.*B(2)-k.*B(1)),(m.*B(2)+k.*B(1).*omega12).*z0.*Mt.*sind(theta2)/(1-Mx.*cosd(theta1)) .* (kx.*CDk/2+ky.*CLk/2).*dz0; %3D matrix, third dimension being the radius
            INT = sum( arg_int , 3); % 2D matrix
            
            
        end
        
        
        A = -rho0.*c0^2.*B2.*sind(theta2)/(8.*pi.*r1/D(rotor_idx).*(1-Mx.*cosd(theta1))); %not sure if it should be theta1 one or theta2 in sin()
        
        
        p = A.*sum(sum(EXP_real.*INT));
        p = A.*sum(sum(sqrt(EXP_real.^2+EXP_imag.^2).*INT));
        
        p_t(t_idx,rotor_idx) = p;
    end
    
    
    
    
end

p_t_total = sum(p_t,2);

p_rms = rms(p_t_total);

dB = 10.*log10(p_rms/p_ref);

end