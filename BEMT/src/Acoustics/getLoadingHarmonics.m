function [CLk,CDk] = getLoadingHarmonics(rotor,dCT,dCP,phi,velocity_dimensional,chord,r,psi,plots,debug)
%{
This is a function for one rotor (doesn't need to be coaxial)
%}

if debug 
    t_res = 30;
else
    t_res = 300;
end

omega = rotor.omega;
R = rotor.R;
Nb = rotor.Nb;

dr = r(1,2)-r(1,1);
dpsi = abs(psi(2,1)-psi(1,1));

rotation_time = 2*pi/omega;

%Could do high resolution here and later do spline interpolation or
%something (when creating the time history).
time_discrete = abs(psi)/omega;

time_interp = repmat(linspace(0,rotation_time,t_res)',1,length(r(1,:)));

%% Transformation from dCT,dCQ/r to dCL,dCD

%see derivation in notebook
%proportionality = rho*A*omega^2*R^2 / (Nb * 0.5*rho*V^2*chord*dr*R*dpsi/(2*pi))
proportionality = 4*pi^2*omega^2*R^3./(Nb*velocity_dimensional.^2.*chord*dr*dpsi);

%inverse rotation matrix (in notebook) -> [cos(phi) sin(phi); -sin(phi) cos(phi)]

%YOU COULD ADD A TURBULENCE FACTOR BY ADDING A RAND MATRIX TO dCT AND dCP
%MATRICES

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

CL_fft = (fft(CL))/length(CL(:,1)); %ABS OR NOT? when I don't use it I at least get a sensible shape for the plot of p(t) for the upper rotor. The reason why goes beyond me.
CD_fft = (fft(CD))/length(CL(:,1));

CL0 = mean(CL);
CLk = [CL0; CL_fft]; %the rest of harmonics

CD0 = mean(CD);
CDk = [CD0; CD_fft]; %the rest of harmonics

CLk = CL_fft;
CDk = CD_fft;

%% Plots

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
    
    %discrete radial stations
    loc = [r(1,1) 0.5 0.9];
    
    
    indx_loc = zeros(1,length(loc));
    str = {};
    i = 1;
    for r_station = loc
        indx = find(abs(r(1,:)-r_station)<dr/2);
        indx_loc(i) = indx;
        str = [str , strcat('r = ' , num2str(r_station))];
        i = i + 1;
    end

    
    %plots with interpolated time
    figure(1)
    plot(time_interp(:,indx_loc),CL(:,indx_loc))
    ylabel('$C_l$ [-]','Interpreter','latex')
    xlabel('$t$ [s]','Interpreter','latex')
    %legend(['r = ', num2str(loc(1))],['r = ', num2str(loc(2))],['r = ', num2str(loc(3))])
    legend(str{:})
    
    figure(2)
    plot(time_interp(:,indx_loc),CD(:,indx_loc))
    ylabel('$C_d$ [-]','Interpreter','latex')
    xlabel('$t$ [s]','Interpreter','latex')
    legend(str{:})
    %plots of harmonics in frequency domain
    
    figure(3)
    plot(0:length(CLk(:,1))-1,CLk(:,indx_loc))
    ylabel('$C_{lk}$ [-]','Interpreter','latex')
    xlabel('Harmonic number','Interpreter','latex')
    legend(str{:})
    
    figure(4)
    plot(0:length(CLk(:,1))-1,CDk(:,indx_loc))
    ylabel('$C_{dk}$ [-]','Interpreter','latex')
    xlabel('Harmonic number','Interpreter','latex')
    legend(str{:})
end

end