function [sound] = hanson_acoustics(coaxial,observer,atm,CLk_struct,CDk_struct,plots)
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


%% User Inputs


harmonics = 20;

approach = 'velocity'; % 'angle' or 'velocity' - they gave a 3dB difference for an example test case. Enough verification!


theta = observer.theta; %[deg] defined as the angle between the rotation axis of the rotor and the line connecting the lower rotor hub with the observer as pictorially depicted in Fig.3
r1 = observer.r1; %[m]
r = observer.r(); %[m]
psi_observer = observer.psi; %[rad]

%% Init

p_ref = 2.*10.^-5; %Pa


% Constants abbreviation - same nomenclature as the paper
rho0 = atm.rho;
c0 = atm.c0;

B(1) = coaxial.rotor(1).Nb;
B(2) = coaxial.rotor(2).Nb;
omega(1) = coaxial.rotor(1).omega;
omega(2)= coaxial.rotor(2).omega;
omega12 = omega(1)./omega(2);
D(1) = 2.*coaxial.rotor(1).R; 
D(2) = 2.*coaxial.rotor(2).R; 
B_D(1) = mean([coaxial.rotor(1).root_chord,coaxial.rotor(1).tip_chord])./D(1);
B_D(2) = mean([coaxial.rotor(2).root_chord,coaxial.rotor(2).tip_chord])./D(2);

BPF(1) = B(1)*omega(1)/(2*pi); %Hz
BPF(2) = B(2)*omega(2)/(2*pi); %Hz


time = linspace(0,1/BPF(1),100); %s


r_vec = linspace(0.01,0.99,20);
psi_vec = linspace(0,2*pi,10);

dpsi = psi_vec(2)-psi_vec(1);
dz0 = r_vec(2)-r_vec(1);

%% Setup matrices

r_stretched = zeros(1,1,length(r_vec));
psi_stretched = zeros(1,1,1,length(psi_vec));

for i = 1:length(r_vec) %can't find a more elegant and quick way of doing this
    r_stretched(1,1,i) = r_vec(i);  
end

%same beun as for r_stretched

for i = 1:length(psi_vec) %can't find a more elegant and quick way of doing this
    psi_stretched(1,1,1,i) = psi_vec(i);  
end

m_range = 0:harmonics-2; %sound harmonic range - START at 0. The -2 is a bit random so that I don't get square matrices
k_range = 0:harmonics-1; %loading harmonic range - START at 0. want to cut off frequencies (unless I apply this Gaussian filter, which for me is the cubic interpolation...)

%dim = [length(k_range),length(m_range),length(r_vec)];

k = repmat(k_range',1,length(m_range),length(r_vec),length(psi_vec));
m = repmat(m_range,length(k_range),1,length(r_vec),length(psi_vec));


z0 = repmat(r_stretched,length(k_range),length(m_range),1,length(psi_vec));
psi = repmat(psi_stretched,length(k_range),length(m_range),length(r_vec),1);


if any(size(k) ~= size(m)) &&  any(size(m) ~= size(z0))  &&  any(size(psi) ~= size(z0))
   error('Sizes of matrices m,k,z0,psi do not match') 
end


%Create a huge matrix

p_k_m_t = zeros(length(k_range),length(m_range),length(time),2); %2 is for the number of rotors


for rotor_idx = [1,2] %1 is upper, 2 is lower rotor
    %% Calculate Clk and Cdk
    
    if rotor_idx == 1
        Clk = CLk_struct.upper;
        Cdk = CDk_struct.upper;
        
    elseif rotor_idx == 2
        Clk = CLk_struct.lower;
        Cdk = CDk_struct.lower;
        
    end
    
    psi = abs(psi);
    if strcmpi(coaxial.rotor(rotor_idx).spin_dir,'CW')
        psi = -1*psi;
    end
    
    %cutoff small harmonics which will not enter the sum. The radius
    %dimension stays the same
    Clk = Clk(1:k_range(end)+1,:);
    Cdk = Cdk(1:k_range(end)+1,:);
    
    for j = 1:length(r_vec)
        CLk(:,:,j) = repmat(Clk(:,j),1,length(m_range));
        CDk(:,:,j) = repmat(Cdk(:,j),1,length(m_range));
    end


    %Create a time-history for each rotor
    for t_idx = 1:length(time)
        t = time(t_idx);
        if B(1) == B(2) && omega(1) == omega(2) % - for a small variation in omega the other equation gave a difference of 2dB
            %% Special case - acoustic interference? In any case, for verification purposes
            
            %S = coaxial.params.interrotor_spacing.*D(rotor_idx)./2;
            
            %match nomenclature
            Omega = omega(1);
            b = B(1);
            
            %% Calculate relative velocities
            
            if strcmpi(approach, 'angle')
                error('Mr not completely fixed with angle method, see how it was done for velocity approach')
                Mx = norm([coaxial.state.axial_vel,coaxial.state.tangent_vel])./c0;
                Mt = Omega.*D(rotor_idx)./2./c0;
                Mr = sqrt(Mx.^2*ones(size(z0))+(Mt*z0).^2);
                
                theta1 = theta-atan2d(coaxial.state.tangent_vel,coaxial.state.axial_vel); %[deg]
                theta2 = theta; % [deg]
                
            elseif strcmpi(approach, 'velocity')
                
                Mx = coaxial.state.axial_vel./c0;
                M_rotation = Omega.*D(rotor_idx)./2./c0;
                MT = coaxial.state.tangent_vel./c0*sin(psi)+M_rotation*ones(size(psi)); %tangential flow component- if I want to add psi here I will have an extra dimension in the integral. TO MAKE THESE CHANGES TO EQ.1 I WOULD LOVE TO UNDERSTAND IT BETTER
                Mt = (Omega.*D(rotor_idx)./2+coaxial.state.tangent_vel)./c0; %tip mach number (due to tangential and rotational flow components)
                
                Mr = sqrt(Mx.^2*ones(size(z0))+MT.^2); %blade section relative mach number
                
                theta1 = theta; %[deg]!!
                theta2 = theta; %[deg]!!
            end
            
            %% Equation 6 exponential (for verification)
            
            arg_exp = (m(:,:,1)-2.*k(:,:,1)).*b.*(psi_observer-pi./2)+m(:,:,1).*b.*(Omega.*r./c0-Omega.*t); %matrix (k_range x m_range)
            EXP_real = -sin(arg_exp); %matrix
            EXP_imag = cos(arg_exp); %matrix
            
            %% Equation 6 integral (for verification)
            
            % I can do element wise on a three-dimensional matrix (k, m, radius)! Maybe I
            % can add time and rotors to the mix - > 5D!
            kx = 2.*b.*Mt./Mr.*(m./(1-Mx.*cosd(theta1))-2.*k).*B_D(rotor_idx);
            ky = -2.*b./(z0.*Mr).*(m.*(Mr.^2.*cosd(theta1)-Mx)./(1-Mx.*cosd(theta1))+2.*k.*Mx).*B_D(rotor_idx); %theta1 is a bit of a guess in the first cosd
            arg_int = Mr.^2.*besselj((m-2.*k).*b,m.*b.*z0.*Mt.*sind(theta2)./(1-Mx.*cosd(theta1))).*(kx.*CDk./2+ky.*CLk./2).*dz0*dpsi; %3D matrix, third dimension being the radius
            INT = sum( arg_int , [3,4]); % 2D matrix
            
        else
            %% Calculate relative velocities
            
            if strcmpi(approach, 'angle')
                Mx = norm([coaxial.state.axial_vel,coaxial.state.tangent_vel])./c0;
                Mt = omega(rotor_idx).*D(rotor_idx)./2./c0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta-atan2d(coaxial.state.tangent_vel,coaxial.state.axial_vel); %[deg]
                theta2 = theta; %[deg]
            elseif strcmpi(approach, 'velocity')
                Mx = coaxial.state.axial_vel./c0;
                Mt = (omega(rotor_idx).*D(rotor_idx)./2+coaxial.state.tangent_vel)./c0;
                Mr = norm([Mx,Mt]);
                
                theta1 = theta; %[deg]
                theta2 = theta; %[deg]
            end
            
            %% Equation 1 exponential (for verification wrt eq.6)
            arg_exp = (m(:,:,1).*B(2)-k(:,:,1).*B(1)).*(psi-pi./2)+(m(:,:,1).*B(2).*omega(2)+k(:,:,1).*B(1).*omega(1)).*(r./c0-t); %matrix (k_range x m_range)
            EXP_real = -sin(arg_exp); %matrix
            EXP_imag = cos(arg_exp); %matrix
            
            
            %% Equation 1 integral (for verification wrt eq.6)
            
            % I can do element wise on a three-dimensional matrix (k, m, radius)! Maybe I
            % can add time and rotors to the mix - > 5D!
            kx = 2.*Mt./Mr.*((m.*B(2)+k.*B(1).*omega12)./(1-Mx.*cosd(theta1))-k.*B(1).*(1+omega12)).*B_D(rotor_idx);
            ky = -2./Mr.*((m.*B(2)+k.*B(1).*omega12).*Mt.^2.*z0.*cosd(theta1)./(1-Mx.*cosd(theta1))-((m.*B(2)-k.*B(1)).*Mx)./z0).*B_D(rotor_idx); %theta1 is a bit of a guess in the first cosd
            arg_int = Mr.^2.*besselj((m.*B(2)-k.*B(1)),(m.*B(2)+k.*B(1).*omega12).*z0.*Mt.*sind(theta2))./(1-Mx.*cosd(theta1)) .* (kx.*CDk./2+ky.*CLk./2).*dz0*dpsi; %3D matrix, third dimension being the radius
            INT = sum( arg_int , [3,4]); % 2D matrix
            
            
        end
        
        
        A = -rho0.*c0.^2.*B(2).*sind(theta2)./(8.*pi.*r1./D(rotor_idx).*(1-Mx.*cosd(theta1))); %not sure if it should be theta1 one or theta2 in sin()
        
        p_k_m_1 = A.*EXP_real(:,:,1,1).*INT;
        p_k_m = A.*sqrt(EXP_real(:,:,1,1).^2+EXP_imag(:,:,1).^2).*INT; % is this what Daniele is referring to? "Yes the imaginary part can be treated as the conventional i*omega*t formulation. In your case write the complex number and then project them in time as you would do with trigonometry"
        
        
        %p = A.*sum(sum(EXP_real(:,:,1).*INT));
        %p = A.*sum(sum(sqrt(EXP_real(:,:,1).^2+EXP_imag(:,:,1).^2).*INT)); 
        
        p_k_m_t(:,:,t_idx,rotor_idx) = p_k_m_1; %k,m,t,rotor
        
    end

end

p_t_upper = sum(p_k_m_t(:,:,:,1),[1,2]); %dimension three is time so we don't want to add in that direction
p_t_lower = sum(p_k_m_t(:,:,:,2),[1,2]); %dimension three is time so we don't want to add in that direction
p_t_total = sum(p_k_m_t,[1,2,4]); %dimension three is time so we don't want to add in that direction

p_rms_upper = rms(p_t_upper);
p_rms_lower = rms(p_t_lower);
p_rms_total = rms(p_t_total);


sound.dB_upper = 20.*log10(p_rms_upper./p_ref);
sound.dB_lower = 20.*log10(p_rms_lower./p_ref);
%sound.dB_total = 20.*log10((p_rms_upper.^2+p_rms_lower.^2)./p_ref.^2); %GIVES AN OVERESTIMATE I think this is what was done in vehicle design tool with thickness and loading noise assuming no acoustic cancellation
sound.dB_total = 20.*log10(p_rms_total./p_ref);

if plots
    %I really have no idea what this is showing me
%     time_mat = repmat(time',1,length(m_range));
%     p_m_t = sum(p_k_m_t,[1,4]);
%     p_m_t = reshape(p_m_t,length(time),length(m_range));
%     figure(10)
%     plot(time_mat(:,1:4),10.*log10(p_m_t(:,1:4)./p_ref))
    
    %% Pressure time history
    figure(11)
    hold on
    pressure_history_total = reshape(sum(p_k_m_t,[1,2,4]),size(time)); %variables already existed
    pressure_history_lower = reshape(sum(p_k_m_t(:,:,:,2),[1,2]),size(time)); %variables already existed
    pressure_history_upper = reshape(sum(p_k_m_t(:,:,:,1),[1,2]),size(time)); %variables already existed
    pressure_history_rms = rms(pressure_history_total); %variable already existed
    plot(time,pressure_history_rms*ones(size(time)))
    plot(time,pressure_history_total)
    plot(time,pressure_history_upper,'-')
    plot(time,pressure_history_lower,'.-')
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$p$ [Pa]','Interpreter','latex')
    legend('p_{RMS}','p(t) total','p(t) upper','p(t) lower')
    
    
    %% dB in frequency spectrum -  trying to replicate figure 5b
    figure(12)
    %dB_spectrum = 20*log10(fft(pressure_history_total)/length(pressure_history_total)/p_ref);
    %plot(dB_spectrum)
    
    p_m_t_total = sum(p_k_m_t,[1,4]);
    p_m_t_total = reshape(p_m_t_total,[length(time),length(m_range)]);
    pressure_harmonics = rms(p_m_t_total);
    freqs = BPF(1)*m_range(2:end);
    semilogx(freqs,20*log10(pressure_harmonics(2:end)/p_ref))
    ylabel('Sound Pressure Level (dB)')
    xlabel('Frequency (Hz)')
    %dB_spectrum = 
end

end