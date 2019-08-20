function [CLk,CDk] = getLoadingHarmonics(rotor,dCT,dCP,phi,velocity_dimensional,chord,r,psi)
%{
This is a function for one rotor (doesn't need to be coaxial)
%}

omega = rotor.omega;
R = rotor.R;
Nb = rotor.Nb;

dr = r(1,2)-r(1,1);
dpsi = abs(psi(2,1)-psi(1,1));

rotation_time = 2*pi/omega;

%Could do high resolution here and later do spline interpolation or
%something (when creating the time history).
time_discrete = abs(psi)/omega;

time_interp = repmat(linspace(0,rotation_time,300)',1,length(r(1,:)));

%% Transformation from dCT,dCQ/r to dCL,dCD

%see derivation in notebook
%proportionality = rho*A*omega^2*R^2 / (Nb * 0.5*rho*V^2*chord*dr*R*dpsi/(2*pi))
proportionality = 4*pi^2*omega^2*R^3./(Nb*velocity_dimensional.^2.*chord*dr*dpsi);

%inverse rotation matrix (in notebook) -> [cos(phi) sin(phi); -sin(phi) cos(phi)]

dCL = proportionality.*(dCT.*cos(phi)+dCP./r.*sin(phi));
dCD = proportionality.*(dCT.*-sin(phi)+dCP./r.*cos(phi));


%% Fourier transform

CL = interp1(time_discrete(:,1),dCL,time_interp(:,1),'PCHIP'); %high res
CD = interp1(time_discrete(:,1),dCD,time_interp(:,1),'PCHIP'); %high res

%this seems to automatically have detected the fundamental frequency. It
%shows as many harmonics as are possible with the sampled points of the
%signal (i.e. the number of rows of CL or the resolution of time_interp).
%What I don't have is CL0 and CD0. I think I can add it by just calculating
%the average. Or.. is the mean already included? It seems to give the same
%result. I don't know if taking abs is correct then

CL_fft = abs(fft(CL))/length(CL(:,1));
CD_fft = abs(fft(CD))/length(CL(:,1));

CL0 = mean(CL);
CLk = [CL0; CL_fft]; %the rest of harmonics

CD0 = mean(CD);
CDk = [CD0; CD_fft]; %the rest of harmonics

CLk = CL_fft;
CDk = CD_fft;

%% Plots
plots = false;

if plots
    %Testing plots - if debug mode on
%     figure(1)
%     plot(rad2deg(abs(psi)),dCL(:,:))
%     ylabel('$C_L$ [-]','Interpreter','latex')
%     xlabel('$\psi$ [deg]','Interpreter','latex')
%     figure(2)
%     plot(rad2deg(abs(psi)),dCD(:,:))
%     ylabel('$C_D$ [-]','Interpreter','latex')
%     xlabel('$\psi$ [deg]','Interpreter','latex')
%     
    %plots wrt time
%     figure(2)
%     plot(time_discrete,dCL(:,:))
%     ylabel('$C_L$ [-]','Interpreter','latex')
%     xlabel('$t$ [s]','Interpreter','latex')
%     figure(2)
%     plot(time_discrete,dCD(:,:))
%     ylabel('$C_D$ [-]','Interpreter','latex')
%     xlabel('$t$ [s]','Interpreter','latex')
%     
    %plots with interpolated time
    figure(1)
    plot(time_interp,CL(:,:))
    ylabel('$C_L$ [-]','Interpreter','latex')
    xlabel('$t$ [s]','Interpreter','latex')
    figure(2)
    plot(time_interp,CD(:,:))
    ylabel('$C_D$ [-]','Interpreter','latex')
    xlabel('$t$ [s]','Interpreter','latex')
    
end

end